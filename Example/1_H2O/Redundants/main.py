from ConcordantModes.options import Options
options_kwargs = {
  'queue'             : "gen3.q,gen4.q,gen6.q,debug.q",
  'program'           : "psi4@master",
  'energyRegex'       : r"Giraffe The Energy is\s+(\-\d+\.\d+)",
  'cartInsert'        : 7,
  'coords'            : "Redundant",
  # 'reducedDisp'       : True,
  #'calc'              : False,
  'successRegex'      : r"beer",
  'off_diag'          : False,
  'off_diag_bands'    : 3,
  'off_diag_limit'    : False,
  'printout_rel_e'    : True,
  'mode_coupling_check' : False 
}
options_obj = Options(**options_kwargs)

# 3. call Mixed Hessian Program
from ConcordantModes.CMA import ConcordantModes
CMA_obj = ConcordantModes(options_obj)
CMA_obj.run()
