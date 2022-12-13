import numpy as np
import os
import re
from numpy.linalg import inv
from numpy import linalg as LA

from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.options import Options
from concordantmodes.s_vectors import SVectors
from concordantmodes.ted import TED
from concordantmodes.zmat import Zmat

np.set_printoptions(precision=9)
os.chdir("./ref_data/f_conv_test/")
FCr = FcRead("fc.dat")
FCr.run()

options = Options()
options.coords = "Redundant"
ZMAT = Zmat(options)
output_test = ZMAT.zmat_read("zmat")
ZMAT.zmat_process(output_test)

ZMAT.zmat_calc()

ZMAT.zmat_compile()

s_vec = SVectors(ZMAT, options, ZMAT.variable_dictionary_init)
s_vec.run(ZMAT.cartesians_init, True)

TED_obj = TED(s_vec.proj, ZMAT)

os.chdir("../../")


def test_f_convert2int():
    os.chdir("./ref_data/f_conv_test/")
    errors = []

    FCint = FcConv(FCr.fc_mat, s_vec, ZMAT, "internal", False, TED_obj, options.units)
    FCint.run()

    FCintR = FcRead("fc_int.dat")
    FCintR.run()

    if np.setdiff1d(FCint.F.round(decimals=10), FCintR.fc_mat).size:
        errors.append(
            "Transformed internal force constants do not match the reference."
        )

    os.chdir("../../")
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


def test_f_convert2cart():
    os.chdir("./ref_data/f_conv_test/")
    errors = []

    FCint = FcConv(FCr.fc_mat, s_vec, ZMAT, "internal", False, TED_obj, options.units)
    FCint.run()
    FCcart = FcConv(FCint.F, s_vec, ZMAT, "cartesian", False, TED_obj, options.units)
    FCcart.run()

    FCintC = FcRead("fc_cart.dat")
    FCintC.run()

    if np.setdiff1d(FCcart.F.round(decimals=10), FCintC.fc_mat).size:
        errors.append(
            "Transformed internal force constants do not match the reference."
        )

    os.chdir("../../")
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


test_f_convert2int()
test_f_convert2cart()
