# 1. build an options object
from ConcordantModes.options import Options
options_kwargs = {
  'queue'            : "gen4.q,gen6.q,gen5.q",
  'program'          : "molpro@2010.1.67+mpi",
  'disp'             : 0.01,
  'energyRegex'      : r"\(T\) energy\s+(\-\d+\.\d+)",
  'cartInsert'       : 11,
  # 'coords'           : "ReDuNdAnT",
  'coords'           : "Custom",
  # 'reducedDisp'       : True,
  # 'geomCheck'        : True,
  # 'calc'             : False,
  'successRegex'     : r"Variable memory released"
}
options_obj = Options(**options_kwargs)

# 2. call Concordant Modes Program
from ConcordantModes.CMA import ConcordantModes
CMA_obj = ConcordantModes(options_obj)
CMA_obj.run()
