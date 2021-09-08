from ConcordantModes.options import Options

options_kwargs = {
    "queue": "gen3.q,gen4.q,gen6.q,debug.q",
    "programInit": "psi4@master",
    "program": "psi4@master",
    "energyRegex": r"Giraffe The Energy is\s+(\-\d+\.\d+)",
    "energyRegexInit": r"Giraffe The Energy is\s+(\-\d+\.\d+)",
    "cartInsertInit": 7,
    "cartInsert": 7,
    "coords": "Redundant",
    # 'reducedDisp'       : True,
    # 'calc'              : False,
    "successRegex": r"beer",
    "successRegexInit": r"beer",
    "off_diag": False,
    "off_diag_bands": 3,
    "off_diag_limit": False,
    "printout_rel_e": False,
    "mode_coupling_check": False,
}
options_obj = Options(**options_kwargs)

# 3. call Mixed Hessian Program
from ConcordantModes.CMA import ConcordantModes

CMA_obj = ConcordantModes(options_obj)
CMA_obj.run()
