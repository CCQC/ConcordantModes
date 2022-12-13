import fileinput
import os
import re
import shutil
import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

from concordantmodes.algorithm import Algorithm
from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.gf_method import GFMethod
from concordantmodes.g_matrix import GMatrix
from concordantmodes.options import Options
from concordantmodes.reap import Reap
from concordantmodes.s_vectors import SVectors
from concordantmodes.ted import TED
from concordantmodes.transf_disp import TransfDisp
from concordantmodes.zmat import Zmat



def test_reap():
    os.chdir("./ref_data/reap_test/")

    options = Options()
    options.cart_insert = 7
    options.coords = "Redundant"

    FC = FcRead("fc.dat")
    FC.run()
    ZMAT = Zmat(options)
    output_test = ZMAT.zmat_read("zmat")
    ZMAT.zmat_process(output_test)
    
    ZMAT.zmat_calc()
    
    ZMAT.zmat_compile()
    s_vec = SVectors(ZMAT, options, ZMAT.variable_dictionary_init)
    s_vec.run(ZMAT.cartesians_init, True)
    
    TED_obj = TED(s_vec.proj, ZMAT)
    f_conv = FcConv(FC.fc_mat, s_vec, ZMAT, "internal", False, TED_obj, options.units)
    f_conv.run()
    
    g_mat = GMatrix(ZMAT, s_vec, options)
    g_mat.run()

    F = np.dot(TED_obj.proj.T, np.dot(f_conv.F, TED_obj.proj))
    G = np.dot(TED_obj.proj.T, np.dot(g_mat.G, TED_obj.proj))
    GF = GFMethod(G, F, options.tol, options.proj_tol, ZMAT, TED_obj)

    GF.run()
    algo = Algorithm(len(GF.L), None, options)
    algo.run()
    
    disps = TransfDisp(
        s_vec,
        ZMAT,
        options.disp,
        GF.L,
        True,
        options.disp_tol,
        TED_obj,
        options,
        algo.indices,
    )
    disps.run()
    
    options.program = "psi4@master"
    prog = options.program
    prog_name = prog.split("@")[0]
    
    options.energy_regex = r"Giraffe The Energy is\s+(\-\d+\.\d+)"
    options.success_regex = r"beer"
    os.chdir("./Disps")
    reap_obj = Reap(
            prog_name,
            ZMAT,
            disps.disp_cart,
            options,
            disps.n_coord,
            GF.L,
            algo.indices,
            options.energy_regex,
            options.gradient_regex,
            options.molly_regex_init,
            options.success_regex
            )
    reap_obj.run()
    
    ref_en = -76.332189646734 
    
    os.chdir("../..")

    assert ref_en == reap_obj.m_en_array[1][1]

test_reap()
