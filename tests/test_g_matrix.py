import os
import shutil
import numpy as np
from concordantmodes.ted import TED
from numpy.linalg import inv
from numpy import linalg as LA

from suite_execute import execute_suite

from concordantmodes.g_matrix import GMatrix
from concordantmodes.options import Options
from concordantmodes.s_vectors import SVectors
from concordantmodes.zmat import Zmat

# s_vec_test = SVectors()


suite = execute_suite("./ref_data/s_vec_test/","Custom",s_vec_bool=True)
suite.run()

# os.chdir("./ref_data/s_vec_test/")
# options = Options()
# options.coords = "Custom"
# ZMAT = Zmat(options)
# output_test = ZMAT.zmat_read("zmat")
# ZMAT.zmat_process(output_test)

# ZMAT.zmat_calc()

# ZMAT.zmat_compile()
# os.chdir("../../")

# s_vec = SVectors(ZMAT, options, ZMAT.variable_dictionary_init)


def test_compute_G():
    errors = []
    suite.s_vec.run(suite.ZMAT.cartesians_init, True)

    g_mat = GMatrix(suite.ZMAT, suite.s_vec, suite.options)
    g_mat.run()

    
    G_ref = [[ 8.00121376e-05, -1.32715756e-05, -1.73338728e-05, -1.73338729e-05, -1.02628585e-05, -1.80496698e-05, -2.12604120e-05, -2.04478488e-05, -2.04478488e-05,  0.00000000e+00,  0.00000000e+00,  5.00626441e-05],
             [-1.32715756e-05,  5.90035575e-04, -1.42964000e-05, -1.42962793e-05,  0.00000000e+00, -1.62927551e-05, -1.62927552e-05,  1.75444777e-05,  1.75444053e-05,  2.59451267e-05, -2.59456007e-05, -3.20145289e-05],
             [-1.73338728e-05, -1.42964000e-05,  5.90035575e-04, -1.47211823e-05,  0.00000000e+00,  7.52213196e-06,  1.84338806e-05, -1.57546251e-05,  2.05416834e-05,  4.34117245e-06,  2.83001423e-05, -1.09112087e-05],
             [-1.73338729e-05, -1.42962793e-05, -1.47211823e-05,  5.90035575e-04,  0.00000000e+00,  7.52323926e-06,  1.84338159e-05,  2.05416834e-05, -1.57546251e-05, -2.82999534e-05, -4.34098373e-06, -1.09111166e-05],
             [-1.02628585e-05,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  5.78617728e-04, -1.21882934e-05, -1.21882934e-05,  5.81936740e-06,  5.82022404e-06,  4.38837729e-06, -4.38818656e-06,  0.00000000e+00],
             [-1.80496698e-05, -1.62927551e-05,  7.52213196e-06,  7.52323926e-06, -1.21882934e-05,  1.91331520e-04,  1.56088852e-05, -7.79557410e-06, -7.79672166e-06, -1.25700962e-05,  1.25695499e-05,  6.80353003e-06],
             [-2.12604120e-05, -1.62927552e-05,  1.84338806e-05,  1.84338159e-05, -1.21882934e-05,  1.56088852e-05,  1.55264109e-04,  1.00850854e-06,  1.00856159e-06, -1.56356767e-05,  1.56366743e-05, -1.42465168e-04],
             [-2.04478488e-05,  1.75444777e-05, -1.57546251e-05,  2.05416834e-05,  5.81936740e-06, -7.79557410e-06,  1.00850854e-06,  1.55208823e-04, -1.12229743e-06, -1.04056524e-05, -2.64304905e-05, -1.16816043e-04],
             [-2.04478488e-05,  1.75444053e-05,  2.05416834e-05, -1.57546251e-05,  5.82022404e-06, -7.79672166e-06,  1.00856159e-06, -1.12229743e-06,  1.55208824e-04,  2.64300382e-05,  1.04052001e-05, -1.16815711e-04],
             [ 0.00000000e+00,  2.59451267e-05,  4.34117245e-06, -2.82999534e-05,  4.38837729e-06, -1.25700962e-05, -1.56356767e-05, -1.04056524e-05,  2.64300382e-05,  3.76136297e-04,  1.94198975e-04, -7.83770679e-05],
             [ 0.00000000e+00, -2.59456007e-05,  2.83001423e-05, -4.34098373e-06, -4.38818656e-06,  1.25695499e-05,  1.56366743e-05, -2.64304905e-05,  1.04052001e-05,  1.94198975e-04,  3.76137316e-04,  7.83768965e-05],
             [ 5.00626441e-05, -3.20145289e-05, -1.09112087e-05, -1.09111166e-05,  0.00000000e+00,  6.80353003e-06, -1.42465168e-04, -1.16816043e-04, -1.16815711e-04, -7.83770679e-05,  7.83768965e-05,  3.78284543e-04]]
    G_ref = np.array(G_ref)
    # print(np.sum(abs(G_ref - g_mat.G)))

    # if np.setdiff1d(np.array(G_ref), np.array(g_mat.G.tolist())).size:
        # errors.append("Computed G-matrix does not match the reference.")

    assert np.sum(abs(G_ref - g_mat.G)) < 1.0e-8


test_compute_G()
