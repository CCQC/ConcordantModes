from concordantmodes.options import Options

options_kwargs = {
    "queue": "gen3.q,gen4.q,gen6.q,debug.q",
    "program_init": "psi4@master",
    "program": "psi4@master",
    "deriv_level_init": 1,
    # "disp" : 0.001,
    "calc_init" : False,
    "energy_regex": r"Giraffe The Energy is\s+(\-\d+\.\d+)",
    "energy_regex_init": [r"Total Gradient",r"tstop"],
    # "energy_regex_init": [r"Total gradient",r"Psi4"],
    "cart_insert_init": 7,
    "cart_insert": 7,
    "coords": "Redundant",
    "success_regex_init": r"beer",
    "success_regex": r"beer",
}
options_obj = Options(**options_kwargs)

# 3. call Concordant Modes Program
from concordantmodes.cma import ConcordantModes

CMA_obj = ConcordantModes(options_obj)
CMA_obj.run()
