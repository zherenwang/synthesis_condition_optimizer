import numpy as np
import re
import copy
import math

__author__ = "Zheren Wang"
__version__ = "1.0"
__maintainer__ = "Zheren Wang"
__email__ = "zherenwang@berkeley.edu"
__status__ = "Development"


class ConditionOptimizer:
    def __init__(self, entries, target_entries, pbx_elts, conc_change=True,
                 V_change=True,
                 pH_change=True):
        self.entries = entries
        self.target_entries = target_entries
        self.pbx_elts = pbx_elts
        self.conc_change = conc_change
        self.V_change = V_change
        self.pH_change = pH_change

    def get_energy_conc(self, entry, pH, V, conc_dict):
        MU_H2O = -2.4583
        PREFAC = 0.0591
        return (
                       entry.uncorrected_energy
                       - MU_H2O * entry.nH2O
                       + entry.npH * PREFAC * pH
                       + entry.nPhi * V
                       + self.get_conc(entry, conc_dict)
               ) * entry.normalization_factor

    def get_conc(self, entry, conc_dict):
        entry_dict = {}
        ion_concentration = []
        entry_list = re.split(r" \+ ", str(entry.name))
        if len(entry_list) == 1:
            entry_dict[entry_list[0]] = [entry.phase_type[0], 1]
        else:
            for i in range(len(entry_list)):
                entry_dict[entry_list[i]] = [entry.phase_type[i], entry.weights[i]]
        for metal in conc_dict:
            for entry in entry_list:
                if re.findall(metal, entry):
                    if entry_dict[entry][0] == "Ion":
                        ion_concentration.append(
                            0.0591 * entry_dict[entry][1] * np.log10(conc_dict[metal])
                        )
        conc_term = sum(ion_concentration)
        return conc_term

    def clean_entries(self, conc_dict):
        target_entries_info = []
        target_entries_ids = []
        for target_entry in self.target_entries:
            target_info = [
                target_entry.npH,
                target_entry.nPhi,
                target_entry.nH2O,
                self.get_conc(target_entry, conc_dict),
            ]
            target_entries_info.append(target_info)
            target_entries_ids.append(target_entry.entry_id)
        for entry in self.entries:
            common_id = [entry_id for entry_id in target_entries_ids if entry_id in entry.entry_id]
            if common_id:
                entry_index = self.entries.index(entry)
                self.entries.pop(entry_index)
                continue

            entry_info = [
                entry.npH,
                entry.nPhi,
                entry.nH2O,
                self.get_conc(entry, conc_dict),
            ]
            for target_info in target_entries_info:
                short_distances = []
                for i in range(4):
                    entry_target_distance = entry_info[i] - target_info[i]
                    if abs(entry_target_distance) < 1e-13:
                        short_distances.append(entry_target_distance)
                if len(short_distances) == 4:
                    entry_index = self.entries.index(entry)
                    self.entries.pop(entry_index)
                    break
        return self.entries


    def energy_convex_hull(self, pH, V, conc_dict):
        all_gs = np.array(
            [self.get_energy_conc(entry, pH, V, conc_dict) for entry in self.entries]
        )
        return np.min(all_gs)

    def target_energy(self, pH, V, conc_dict):
        target_energies = [
            self.get_energy_conc(entry, pH, V, conc_dict)
            for entry in self.target_entries
        ]
        return np.min(target_energies)

    def energy_above_hull(self, pH, V, conc_dict):
        return self.target_energy(pH, V, conc_dict) - self.energy_convex_hull(
            pH, V, conc_dict
        )

    def energy_target_below_others(self, pH, V, conc_dict):
        return self.energy_convex_hull(pH, V, conc_dict) - self.target_energy(
            pH, V, conc_dict
        )

    def get_gradient(self, pH, V, conc_dict):
        old_conc_dict = copy.deepcopy(conc_dict)
        self.gradient = {}
        self.gradient["gradient_V"] = (
                                              self.energy_target_below_others(pH, V + 0.01, conc_dict)
                                              - self.energy_target_below_others(pH, V - 0.01, conc_dict)
                                      ) / (2 * 0.01)
        self.gradient["gradient_pH"] = (
                                               self.energy_target_below_others(pH + 0.01, V, conc_dict)
                                               - self.energy_target_below_others(pH - 0.01, V, conc_dict)
                                       ) / (2 * 0.01)
        for element in conc_dict:
            name = "gradient_" + str(element)
            conc_dict[element] += 0.01
            self.gradient[name] = (
                                          self.energy_target_below_others(pH, V, conc_dict)
                                          - self.energy_target_below_others(pH, V, old_conc_dict)
                                  ) / 0.01
            conc_dict[element] = old_conc_dict[element]

        return self.gradient

    def get_thermodynamic_competition(self, conc_dict, V, pH):
        return -1 * self.energy_target_below_others(pH, V, conc_dict)

    def optimizer(self, conc_dict, V, pH, conc_limit, iter_max, verbose=False):
        # RMSProp
        learning_rate = 0.01
        gamma = 0.99
        pH_list = [V]
        V_list = [pH]
        entry_and_energy = {}
        iter_count = 0
        eps = 1e-6
        self.gradient = self.get_gradient(pH, V, conc_dict)
        step = {"V": 0, "pH": 0}
        for element in conc_dict.keys():
            step[str(element)] = 0
        while any(self.gradient.values()) and iter_count < iter_max:
            self.gradient = self.get_gradient(pH, V, conc_dict)
            if verbose:
                print(self.gradient)

            if self.gradient["gradient_V"] != 0 and self.V_change and V <= 2 and V >= -2:
                step["V"] = (
                        gamma * step["V"] + (1 - gamma) * self.gradient["gradient_V"] ** 2
                )
                V += 0.1 * (
                        learning_rate
                        / math.sqrt(step["V"] + eps)
                        * self.gradient["gradient_V"]
                )

            if self.gradient["gradient_pH"] != 0 and self.pH_change and pH <= 14 and pH >= 0:
                step["pH"] = (
                        gamma * step["pH"] + (1 - gamma) * self.gradient["gradient_pH"] ** 2
                )
                pH += (
                        learning_rate
                        / math.sqrt(step["pH"] + eps)
                        * self.gradient["gradient_pH"]
                )
            if self.conc_change:
                for element in conc_dict:
                    name = "gradient_" + str(element)
                    if self.gradient[name] != 0 and conc_dict[element] < conc_limit:
                        step[str(element)] = (
                                gamma * step[str(element)]
                                + (1 - gamma) * learning_rate * self.gradient[name] ** 2
                        )
                        conc_dict[element] += 0.06 * (
                                learning_rate
                                / math.sqrt(step[str(element)] + eps)
                                * self.gradient[name]
                        )
            all_gs = np.array(
                [self.get_energy_conc(entry, pH, V, conc_dict) for entry in self.entries]
            )
            if verbose:
                print(self.entries[np.argmin(all_gs)])
            entry_and_energy = dict(zip([entry.name for entry in self.entries], all_gs))
            entry_and_energy[self.target_entries[0].name] = self.target_energy(pH, V, conc_dict)
            if verbose:
                print('thermodynamic competition is ' + str(-1 * self.energy_target_below_others(pH, V, conc_dict)))
            pH_list.append(pH)
            V_list.append(V)
            iter_count += 1
            if verbose:
                print(iter_count)
        if verbose:
            print('pH = ' + str(pH))
            print('V =' + str(V) + ' V')
            print('conc_dict is ', conc_dict)
            print('thermodynamic competition is ' + str(-1 * self.energy_target_below_others(pH, V, conc_dict)))
        return pH, V, conc_dict, pH_list, V_list, entry_and_energy, -1 * self.energy_target_below_others(pH, V,
                                                                                                         conc_dict)
