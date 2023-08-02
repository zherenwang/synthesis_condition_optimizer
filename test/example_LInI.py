from synthesis_condition_optimizer.synthesis_condition_optimizer import ConditionOptimizer
from pymatgen.core.composition import Composition
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
import datetime


def get_target(target, entries):
    entry_dict = {}
    for entry in entries:
        if target in entry.name:
            entry_dict[entry] = entry.normalized_energy
    if entry_dict:
        sort_entry = sorted(entry_dict.items(), key=lambda item: item[1])
        return sort_entry[0]


def get_V(pH):
    return 1.626 - 0.0591 * pH


if __name__ == '__main__':

    start_time = datetime.datetime.now()

    target_id = []

    target = "LiIn(IO3)4"
    conc_dict = {'Li': 1.5, 'In': 0.1, 'I': 0.8}
    pbx_elts = list(conc_dict.keys())
    pH = 0.66
    V = get_V(pH)
    print(V)
    print('target is ' + target)
    comp = Composition(target).as_dict()
    comp_dict = {}
    for elt in pbx_elts:
        comp_dict[elt] = comp[elt]

    print("Getting Entries")

    MPR = MPRester("Your API key")
    entries = MPR.get_pourbaix_entries(pbx_elts)

    pd = PourbaixDiagram(entries, conc_dict=conc_dict, comp_dict=comp_dict, filter_solids=True)
    all_entries = pd._processed_entries
    print(len(pd._unprocessed_entries))
    target_entry = get_target(target, entries)

    if target_entry:
        target_entry = target_entry[0]
        target_entries = [target_entry]
        target_id = target_entries[0].entry_id
        cvx = ConditionOptimizer(all_entries, target_entries, pbx_elts, conc_change=False, V_change=False,
                                 pH_change=False)
        cvx.clean_entries(conc_dict)
        # pH, V, conc_dict, pH_list, V_list, entry_and_energy, thermodynamic_competition = cvx.optimizer(
        #     conc_dict, V=get_V(pH), pH=pH, conc_limit=1.5, iter_max=1
        # )
        thermodynamic_competition = cvx.get_thermodynamic_competition(conc_dict, V=get_V(pH), pH=pH)

        print('pH = ' + str(pH))
        print('redox potential = ' + str(V) + ' V')
        print('conc_dict is ', conc_dict, '(Unit: mol / L)')
        print('thermodynamic competition is ' + str(thermodynamic_competition) + ' eV/atom')

        end_time = datetime.datetime.now()
        time_interval = (end_time - start_time).seconds
        print('running time is ' + str(time_interval) + ' second')

        print("Done!")
