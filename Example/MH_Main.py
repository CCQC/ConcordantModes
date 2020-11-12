from MixedHessian.options import Options
options_kwargs = {
  'queue'             : "gen3.q,gen4.q,gen6.q,debug.q",
  'program'           : "psi4@master",
  'energyRegex'       : r"Giraffe The Energy is\s+(\-\d+\.\d+)",
  'cartInsert'        : 7,
  'successRegex'      : r"beer"
}
options_obj = Options(**options_kwargs)

# 3. call Mixed Hessian Program
from MixedHessian.MH import MixedHessian
MH_obj = MixedHessian(options_obj)
MH_obj.run()
