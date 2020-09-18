# 1. build a submit function
# import subprocess as sp
# import shutil
# import os
# 0. Copy input hessian
# shutil.copy('output.default.hess','output.default.{}.hess'.format(os.getpid()))

# 1. Define submit function
# def submit(optns):
    # sp.run(["python2", "/home/vulcan/mel64643/bin/optavc/submitter.py", optns.queue, optns.program, str(optns.job_array_range[1]), str(optns.resub)])

# 2. build an options object
from MixedHessian.options import Options
options_kwargs = {
  # 'input_name'        : "input.dat",
  # 'output_name'       : "output.dat",
  'queue'             : "gen3.q,gen4.q,gen6.q,debug.q",
  'program'           : "psi4@master",
  'adisp'             : 0.005,
  'basis'             : "cc-pVTZ"
  # 'program'           : "molpro@2010.1.67+mpi",
  # 'template_file_path': "template.dat",
  # 'energy_regex'      : r"\!RHF\-UCCSD\(T\) energy\s+(\-\d+\.\d+)",
  # 'energy_regex'      : r"\(T\) energy\s+(\-\d+\.\d+)",
  # 'success_regex'     : r"Variable memory released" ,
  # 'submitter'         : submit,
  # 'job_title'         : 'min3 optimization',
  # 'email_from'        : 'mel64643@uga.edu',
  # 'email_to'          : 'mitchlahm@gmail.com',
  # 'job_array'         : True,
  # 'resub'             : True,
  # 'maxiter'           : 40,
}
options_obj = Options(**options_kwargs)

# 3. call Mixed Hessian Program
from MixedHessian.MH import MixedHessian
MH_obj = MixedHessian(options_obj)
MH_obj.run()
