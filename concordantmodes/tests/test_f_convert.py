import numpy as np
import os
import re
from numpy.linalg import inv
from numpy import linalg as LA

from suite_execute import execute_suite

from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.options import Options
from concordantmodes.s_vectors import SVectors
from concordantmodes.ted import TED
from concordantmodes.zmat import Zmat

np.set_printoptions(precision=9)

suite = execute_suite("./ref_data/f_conv_test/", "Redundant")
suite.run()


def test_f_convert2int():
    errors = []

    FCint = FcConv(
        suite.FC.fc_mat,
        suite.s_vec,
        suite.ZMAT,
        "internal",
        False,
        suite.TED_obj,
        suite.options.units,
    )
    FCint.run()

    FCintR = FcRead(suite.path + "/fc_int.dat")
    FCintR.run()

    if np.setdiff1d(FCint.F.round(decimals=10), FCintR.fc_mat).size:
        errors.append(
            "Transformed internal force constants do not match the reference."
        )

    assert not errors, "errors occured:\n{}".format("\n".join(errors))


def test_f_convert2cart():
    errors = []

    FCint = FcConv(
        suite.FC.fc_mat,
        suite.s_vec,
        suite.ZMAT,
        "internal",
        False,
        suite.TED_obj,
        suite.options.units,
    )
    FCint.run()
    FCcart = FcConv(
        FCint.F,
        suite.s_vec,
        suite.ZMAT,
        "cartesian",
        False,
        suite.TED_obj,
        suite.options.units,
    )
    FCcart.run()

    FCintC = FcRead(suite.path + "/fc_cart.dat")
    FCintC.run()

    if np.setdiff1d(FCcart.F.round(decimals=10), FCintC.fc_mat).size:
        errors.append(
            "Transformed internal force constants do not match the reference."
        )

    assert not errors, "errors occured:\n{}".format("\n".join(errors))


# test_f_convert2int()
# test_f_convert2cart()
