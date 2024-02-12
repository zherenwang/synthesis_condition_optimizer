# synthesis_condition_optimizer
This repository includes samples codes that calculate the thermodynamic competition and optimize synthesis conditions by minimizing thermodynamic competition.
Version=1.0

installation of the synthesis condition optimizer (only requires pre-install numpy) and pymatgen
~~~
pip install -e.
pip install pymatgen
~~~


Examples:
Calculate the thermodynamic competition for LiIn(IO3)4
~~~
An example is in test/example_LInI.py
~~~
Output (test on apple M2 chip)
~~~
pH = 0.66
redox potential = 1.5869939999999998 V
conc_dict is  {'Li': 1.5, 'In': 0.1, 'I': 0.8} (Unit: mol / L)
thermodynamic competition is -0.056 eV/atom
running time is 1 second
Done!
~~~

Calculate the thermodynamic competition for LiFePO4
~~~
An example is in test/example_LFP.py
~~~
Output (test on apple M2 chip)
~~~
pH = 8.29
redox potential = -0.506239 V
conc_dict is  {'Li': 0.75, 'Fe': 0.25, 'P': 0.28} (Unit: mol / L)
thermodynamic competition is -0.058 eV/atom
running time is 5 second
Done!
~~~

Optimize synthesis condition to synthesis BaCO3

~~~
An example is in test/example_optimization.py
~~~
Output (test on apple M2 chip)
~~~
pH = 9.9
redox potential = 0 V
conc_dict is  {'Ba': 2.00, 'C': 2.00} (Unit: mol / L)
thermodynamic competition is -0.38 eV/atom
running time is 11 second
Done!

~~~

If you find the codes and data useful, please consider citing our paper:
~~~

@article{wang_2024_MTC,
	author = {Wang, Zheren and Sun, Yingzhi and Cruse, Kevin and Zeng, Yan and Fei, Yuxing and Liu, Zexuan and Shangguan, Junyi and Byeon, Young-Woon and Jun, KyuJung and He, Tanjin and Sun, Wenhao and Ceder, Gerbrand},
	doi = {10.1038/s44160-023-00479-0},
	journal = {Nature Synthesis},
	title = {Optimal thermodynamic conditions to minimize kinetic by-products in aqueous materials synthesis},
	year = {2024}
	}
~~~