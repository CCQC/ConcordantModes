import os
import shutil
import numpy as np
from concordantmodes.ted import TED
from numpy.linalg import inv
from numpy import linalg as LA

from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.gf_method import GFMethod
from concordantmodes.g_matrix import GMatrix
from concordantmodes.options import Options
from concordantmodes.s_vectors import SVectors
from concordantmodes.ted import TED
from concordantmodes.zmat import Zmat


os.chdir("./ref_data/f_read_test/")
options = Options()
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

os.chdir("../../")


def test_gf_method():
    errors = []
    F = np.dot(TED_obj.proj.T, np.dot(f_conv.F, TED_obj.proj))
    G = np.dot(TED_obj.proj.T, np.dot(g_mat.G, TED_obj.proj))
    GF = GFMethod(G, F, options.tol, options.proj_tol, ZMAT, TED_obj)

    GF.run()

    ref_freq = [
        347.2028653723603,
        1073.4528968170207,
        1099.7617262833578,
        1171.279871709732,
        1405.105997416581,
        1483.2609607580182,
        1489.2922133701677,
        1507.312943326969,
        2998.885895844555,
        3054.76300252422,
        3134.79345416958,
        3837.8500878962245,
    ]

    if np.setdiff1d(
        GF.freq.round(decimals=10), np.array(ref_freq).round(decimals=10)
    ).size:
        errors.append(
            "Frequencies computed via the GF method do not match the reference."
        )

    assert not errors, "errors occured:\n{}".format("\n".join(errors))


test_gf_method()
