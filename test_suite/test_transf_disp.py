import os
import shutil
import numpy as np
from concordantmodes.ted import TED
from numpy.linalg import inv
from numpy import linalg as LA

from concordantmodes.algorithm import Algorithm
from concordantmodes.f_convert import FcConv
from concordantmodes.f_read import FcRead
from concordantmodes.gf_method import GFMethod
from concordantmodes.g_matrix import GMatrix
from concordantmodes.options import Options
from concordantmodes.s_vectors import SVectors
from concordantmodes.ted import TED
from concordantmodes.transf_disp import TransfDisp
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

F = np.dot(TED_obj.proj.T, np.dot(f_conv.F, TED_obj.proj))
G = np.dot(TED_obj.proj.T, np.dot(g_mat.G, TED_obj.proj))
GF = GFMethod(G, F, options.tol, options.proj_tol, ZMAT, TED_obj)

GF.run()
algo = Algorithm(len(GF.L), None, options)
algo.run()

os.chdir("../../")


def test_transf_disp():
    errors = []

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

    disp_ref = [
        [-1.3750734515947842, -0.024279474843673623, 0.0029358198799196837],
        [1.3060102262375441, 0.1207018145868287, -0.0023356196755552655],
        [-2.0778682947484093, 1.9096086118030944, -0.0008923821222545903],
        [-2.1129906861119383, -0.9795514202073832, 1.683039582462003],
        [-2.1049748848658, -0.979121064823173, -1.6808319627150445],
        [1.941172931083388, -1.5774905865156932, -0.0019238478290682895],
    ]

    # print(disps.p_disp[3][3].tolist())

    if np.setdiff1d(np.array(disp_ref), disps.p_disp[3][3].tolist()).size:
        errors.append("Computed Displacement does not match the reference.")

    assert not errors, "errors occured:\n{}".format("\n".join(errors))


test_transf_disp()
