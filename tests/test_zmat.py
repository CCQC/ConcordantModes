import pytest
import numpy as np
import os
import shutil
import re

from concordantmodes.int2cart import Int2Cart
from concordantmodes.transf_disp import TransfDisp
from concordantmodes.options import Options
from concordantmodes.zmat import Zmat

coord1 = "Redundant"
coord2 = "ZMAT"
coord3 = "Custom"
file1 = "zmat_red"
file2 = "zmat_zmat"
file3 = "zmat_custom"
# ZMAT read data
ref_Red = ["1 2\n", "1 3\n", "1 4\n", "1 5\n", "2 6\n"]
ref_ZMAT = ["C\n", "O 1\n", "H 1 2\n", "H 1 2 3\n", "H 1 2 4\n", "H 2 1 3\n"]
ref_Custom = [
        "1 2\n",
        "1 3\n",
        "1 4\n",
        "1 5\n",
        "2 6\n",
        "2 1 3\n",
        "2 1 4\n",
        "2 1 5\n",
        "3 1 4\n",
        "4 1 5\n",
        "5 1 3\n",
        "6 2 1\n",
        "6 2 1 4 T\n",
        "4 1 3 5 O\n",
        "3 1 2 6 Lx\n",
        "3 1 2 6 Ly\n",
        "5 1 2 4 L\n",
    ]
zmat_read = [(coord1,ref_Red,file1),(coord2,ref_ZMAT,file2),(coord2,ref_Custom,file3)]

@pytest.mark.parametrize(
    "option, expected, file_name",zmat_read 
)

def test_zmat_read(option, expected, file_name):
    os.chdir("./ref_data/zmat_test/")
    options = Options()

    options.coords = option
    ZMAT = Zmat(options)

    output_test = ZMAT.zmat_read(file_name)

    os.chdir("../../")
    assert expected == output_test


def test_zmat_process():
    os.chdir("./ref_data/zmat_test/")
    options = Options()
    errors = []

    options.coords = "ZMAT"
    ZMAT = Zmat(options)
    output_test_zmat = ZMAT.zmat_read("zmat_zmat")
    ZMAT.zmat_process(output_test_zmat)
    ref_bond_indices = [["2", "1"], ["3", "1"], ["4", "1"], ["5", "1"], ["6", "2"]]
    ref_angle_indices = [
        ["3", "1", "2"],
        ["4", "1", "2"],
        ["5", "1", "2"],
        ["6", "2", "1"],
    ]
    ref_torsion_indices = [
        ["4", "1", "2", "3"],
        ["5", "1", "2", "4"],
        ["6", "2", "1", "3"],
    ]
    ref_bond_variables = ["R1", "R2", "R3", "R4", "R5"]
    ref_angle_variables = ["A2", "A3", "A4", "A5"]
    ref_torsion_variables = ["D3", "D4", "D5"]
    if not ref_bond_indices == ZMAT.bond_indices:
        errors.append("ZMAT bond indices")
    if not ref_angle_indices == ZMAT.angle_indices:
        errors.append("ZMAT angle indices")
    if not ref_torsion_indices == ZMAT.torsion_indices:
        errors.append("ZMAT torsion indices")
    if not ref_bond_variables == ZMAT.bond_variables:
        errors.append("ZMAT bond variables")
    if not ref_angle_variables == ZMAT.angle_variables:
        errors.append("ZMAT angle variables")
    if not ref_torsion_variables == ZMAT.torsion_variables:
        errors.append("ZMAT torsion variables")

    options.coords = "Redundant"
    ZMAT = Zmat(options)
    output_test_red = ZMAT.zmat_read("zmat_red")
    ZMAT.zmat_process(output_test_red)
    ref_bond_indices = np.array(
        [["1", "2"], ["1", "3"], ["1", "4"], ["1", "5"], ["2", "6"]]
    )
    ref_angle_indices = np.array(
        [
            ["2", "1", "3"],
            ["2", "1", "4"],
            ["2", "1", "5"],
            ["1", "2", "6"],
            ["3", "1", "4"],
            ["3", "1", "5"],
            ["4", "1", "5"],
        ]
    )
    ref_torsion_indices = np.array(
        [
            ["3", "1", "2", "4"],
            ["3", "1", "2", "5"],
            ["3", "1", "2", "6"],
            ["2", "1", "3", "4"],
            ["2", "1", "3", "5"],
            ["4", "1", "2", "5"],
            ["4", "1", "2", "6"],
            ["2", "1", "4", "3"],
            ["2", "1", "4", "5"],
            ["5", "1", "2", "6"],
            ["2", "1", "5", "3"],
            ["2", "1", "5", "4"],
            ["4", "1", "3", "5"],
            ["3", "1", "4", "5"],
            ["3", "1", "5", "4"],
        ]
    )
    ref_bond_variables = ["R1", "R2", "R3", "R4", "R5"]
    ref_angle_variables = ["A1", "A2", "A3", "A4", "A5", "A6", "A7"]
    ref_torsion_variables = [
        "D1",
        "D2",
        "D3",
        "D4",
        "D5",
        "D6",
        "D7",
        "D8",
        "D9",
        "D10",
        "D11",
        "D12",
        "D13",
        "D14",
        "D15",
    ]
    # print(np.setdiff1d(ref_bond_indices,ZMAT.bond_indices))
    if np.setdiff1d(ref_bond_indices, ZMAT.bond_indices).size:
        errors.append("Redundant bond indices")
    if np.setdiff1d(ref_angle_indices, ZMAT.angle_indices).size:
        errors.append("Redundant angle indices")
    if np.setdiff1d(ref_torsion_indices, ZMAT.torsion_indices).size:
        errors.append("Redundant torsion indices")
    if np.setdiff1d(ref_bond_variables, ZMAT.bond_variables).size:
        errors.append("Redundant bond variables")
    if np.setdiff1d(ref_angle_variables, ZMAT.angle_variables).size:
        errors.append("Redundant angle variables")
    if np.setdiff1d(ref_torsion_variables, ZMAT.torsion_variables).size:
        errors.append("Redundant torsion variables")

    options.coords = "Custom"
    ZMAT = Zmat(options)
    output_test_red = ZMAT.zmat_read("zmat_custom")
    ZMAT.zmat_process(output_test_red)
    ref_bond_indices = [("1", "2"), ("1", "3"), ("1", "4"), ("1", "5"), ("2", "6")]
    ref_angle_indices = [
        ("2", "1", "3"),
        ("2", "1", "4"),
        ("2", "1", "5"),
        ("3", "1", "4"),
        ("4", "1", "5"),
        ("5", "1", "3"),
        ("6", "2", "1"),
    ]
    ref_torsion_indices = [("6", "2", "1", "4")]
    ref_oop_indices = [("4", "1", "3", "5")]
    ref_lin_indices = [("5", "1", "2", "4")]
    ref_linx_indices = [("3", "1", "2", "6")]
    ref_liny_indices = [("3", "1", "2", "6")]
    ref_bond_variables = ["R1", "R2", "R3", "R4", "R5"]
    ref_angle_variables = ["A1", "A2", "A3", "A4", "A5", "A6", "A7"]
    ref_torsion_variables = ["D1"]
    ref_oop_variables = ["O1"]
    ref_lin_variables = ["L1"]
    ref_linx_variables = ["Lx1"]
    ref_liny_variables = ["Ly1"]

    if np.setdiff1d(ref_bond_indices, ZMAT.bond_indices).size:
        errors.append("Redundant bond indices")
    if np.setdiff1d(ref_angle_indices, ZMAT.angle_indices).size:
        errors.append("Redundant angle indices")
    if np.setdiff1d(ref_torsion_indices, ZMAT.torsion_indices).size:
        errors.append("Redundant torsion indices")
    if np.setdiff1d(ref_oop_indices, ZMAT.oop_indices).size:
        errors.append("Redundant out-of-plane indices")
    if np.setdiff1d(ref_lin_indices, ZMAT.lin_indices).size:
        errors.append("Redundant lin indices")
    if np.setdiff1d(ref_linx_indices, ZMAT.linx_indices).size:
        errors.append("Redundant linx indices")
    if np.setdiff1d(ref_liny_indices, ZMAT.liny_indices).size:
        errors.append("Redundant liny indices")
    if not ref_bond_variables == ZMAT.bond_variables:
        errors.append("Redundant bond variables")
    if not ref_angle_variables == ZMAT.angle_variables:
        errors.append("Redundant angle variables")
    if not ref_torsion_variables == ZMAT.torsion_variables:
        errors.append("Redundant torsion variables")
    if not ref_oop_variables == ZMAT.oop_variables:
        errors.append("Redundant out-of-plane variables")
    if not ref_lin_variables == ZMAT.lin_variables:
        errors.append("Redundant lin variables")
    if not ref_linx_variables == ZMAT.linx_variables:
        errors.append("Redundant linx variables")
    if not ref_linx_variables == ZMAT.linx_variables:
        errors.append("Redundant linx variables")

    os.chdir("../../")

    assert not errors, "errors occured:\n{}".format("\n".join(errors))


# Only need to test the Custom internal coordinates
def test_zmat_calc():
    os.chdir("./ref_data/zmat_test/")
    options = Options()
    errors = []

    options.coords = "Custom"
    ZMAT = Zmat(options)
    output_test_red = ZMAT.zmat_read("zmat_custom")
    ZMAT.zmat_process(output_test_red)

    ZMAT.zmat_calc()

    var_dict_ref = {
        "R1": 2.685006407296404,
        "R2": 2.0576342503795035,
        "R3": 2.068739344106188,
        "R4": 2.068739340043222,
        "R5": 1.8130883453288267,
        "A1": 106.8765896566676,
        "A2": 112.2824413394013,
        "A3": 112.28244150454864,
        "A4": 108.22376477570377,
        "A5": 108.7851835780631,
        "A6": 108.22360552816234,
        "A7": 107.41157765289944,
        "D1": 61.480653653479756,
        "O1": 57.2185581514689,
        "L1": -45.926274971749756,
        "Lx1": -54.67048011096265,
        "Ly1": -0.0022793163188044303,
    }
    var_dict_custom = ZMAT.variable_dictionary_final

    if np.setdiff1d(var_dict_ref, var_dict_custom).size:
        errors.append("Custom variables do not match.")

    os.chdir("../../")
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


# Only need to test the Custom internal coordinates
def test_zmat_compile():
    os.chdir("./ref_data/zmat_test/")
    options = Options()
    errors = []

    options.coords = "Custom"
    ZMAT = Zmat(options)
    output_test_custom = ZMAT.zmat_read("zmat_custom")
    ZMAT.zmat_process(output_test_custom)

    ZMAT.zmat_calc()

    ZMAT.zmat_compile()

    print(ZMAT.index_dictionary)

    index_dict_ref = {
        "R1": ("1", "2"),
        "R2": ("1", "3"),
        "R3": ("1", "4"),
        "R4": ("1", "5"),
        "R5": ("2", "6"),
        "A1": ("2", "1", "3"),
        "A2": ("2", "1", "4"),
        "A3": ("2", "1", "5"),
        "A4": ("3", "1", "4"),
        "A5": ("4", "1", "5"),
        "A6": ("5", "1", "3"),
        "A7": ("6", "2", "1"),
        "D1": ("6", "2", "1", "4"),
        "O1": ("4", "1", "3", "5"),
        "L1": ("5", "1", "2", "4"),
        "Lx1": ("3", "1", "2", "6"),
        "Ly1": ("3", "1", "2", "6"),
    }
    index_dict_custom = ZMAT.index_dictionary

    if np.setdiff1d(index_dict_ref, index_dict_custom).size:
        errors.append("Custom indices do not match.")

    os.chdir("../../")
    assert not errors, "errors occured:\n{}".format("\n".join(errors))


# test_zmat_read_ZMAT()
