class Options(object):
    def __init__(self, **kwargs):
        self.benchmark_full = kwargs.pop("benchmark_full", False)
        self.bond_tol = kwargs.pop("inter_tol", 3.778)
        self.calc = kwargs.pop("calc", True)
        self.cart_coords = kwargs.pop("cart_coords", "Bohr")
        self.calc_init = kwargs.pop("calc_init", True)
        self.cart_insert = kwargs.pop("cart_insert", -1)
        self.cart_insert_init = kwargs.pop("cart_insert_init", -1)
        self.clean_house = kwargs.pop("clean_house", True)
        self.cluster = kwargs.pop("cluster", "vulcan")
        self.coords = kwargs.pop("coords", "ZMAT")
        self.deriv_level = kwargs.pop("deriv_level", 0)
        self.deriv_level_init = kwargs.pop("deriv_level_init", 0)
        self.dir_reap = kwargs.pop("dir_reap", True)
        self.disp = kwargs.pop("disp", 0.01)
        self.disp_check = kwargs.pop("disp_check", False)
        self.disp_tol = kwargs.pop("disp_tol", 0.0001)
        self.energy_regex = kwargs.pop("energy_regex", "")
        self.energy_regex_init = kwargs.pop("energy_regex_init", "")
        self.gen_disps = kwargs.pop("gen_disps",True)
        self.gen_disps_init = kwargs.pop("gen_disps_init",True)
        self.geom_check = kwargs.pop("geom_check", False)
        self.gradient_regex = kwargs.pop("gradient_regex", "")
        self.interatomic_distance = kwargs.pop("interatomic_distance", False)
        self.man_proj = kwargs.pop("man_proj", False)
        self.mode_coupling_check = kwargs.pop("mode_coupling_check", False)
        self.molly_regex_init = kwargs.pop("molly_regex_init", "")
        self.nslots = kwargs.pop("nslots", 1)
        self.off_diag = kwargs.pop("off_diag", False)
        self.off_diag_bands = kwargs.pop("off_diag_bands", 1)
        self.off_diag_limit = kwargs.pop("off_diag_limit", False)
        self.printout_rel_e = kwargs.pop("printout_rel_e", True)
        self.program = kwargs.pop("program", "molpro@2010.1.67+mpi")
        self.program_init = kwargs.pop("program_init", "molpro@2010.1.67+mpi")
        self.proj_tol = kwargs.pop("proj_tol", 1.0e-14)
        self.queue = kwargs.pop("queue", "gen4.q")
        self.reduced_disp = kwargs.pop("reduced_disp", False)
        self.second_order = kwargs.pop("second_order", False)
        self.rmsd = kwargs.pop("rmsd", False)
        self.spin = kwargs.pop("spin", 1)
        self.success_regex = kwargs.pop("success_regex", "")
        self.success_regex_init = kwargs.pop("success_regex_init", "")
        self.tight_disp = kwargs.pop("tight_disp", False)
        self.tol = kwargs.pop("tol", 1.0e-14)
        self.units = kwargs.pop("units", "HartreeBohr")
