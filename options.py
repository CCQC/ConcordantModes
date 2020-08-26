class Options(object):
    def __init__(self, **kwargs):
        self.queue = kwargs.pop("queue", "gen4.q")
        self.nslots = kwargs.pop("nslots", 4)
        self.program = kwargs.pop('program', "molpro@2010.1.67+mpi")
        self.rdisp = kwargs.pop('rdisp', 0.005)
        self.adisp = kwargs.pop('adisp', 0.01)
        self.basis = kwargs.pop('basis', "cc-pVDZ")

        ## These options may be helpful in the future, or maybe not!
        # self.input_name = kwargs.pop("input_name", "input.dat")
        # self.output_name = kwargs.pop("output_name", "output.dat")
        # self.job_array_range = kwargs.pop("job_array", False)  # needs to be set by calling function
        # self.template_file_path = kwargs.pop("template_file_path", "template.dat")
        # self.files_to_copy = kwargs.pop("files_to_copy", [])
        # self.input_units = kwargs.pop("input_units", "bohr")
        # self.point_group = kwargs.pop("point_group", None)
        # self.maxiter = kwargs.pop("maxiter", 20)
        # self.job_array = kwargs.pop("job_array", False)
        # self.email = kwargs.pop("email", None)
        # self.email_opts = kwargs.pop("email_opts", 'ae')
        # self.memory = kwargs.pop("memory", "")
        # self.name = kwargs.pop("name", "")
        # self.resub = kwargs.pop("resub",False)
        # self.resub_test = kwargs.pop("resub_test",False)
        
        ## These options may be helpful in the future for porting over to Sapelo
        # self.cluster = kwargs.pop("cluster", "").upper()
        # self.wait_time = kwargs.pop("wait_time", None)
        # self.time_limit = kwargs.pop("time_limit", None)
