from . import test_adaptive as tk
import pele_platform.main as main
import pele_platform.Errors.custom_errors as ce
import pele_platform.constants.constants as cs
import glob
import os


test_path = os.path.join(cs.DIR, "Examples")
METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_metals.yaml")
NO_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_no_metal_constraints.yaml")
FAIL_PERMISSIVE_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/fail_input_permissive_constraints.yaml")
PASS_PERMISSIVE_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/pass_input_permissive_constraints.yaml")
ALL_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_all_metal_constraints.yaml")
PERMISSIVE_EXCEPTION = os.path.join(test_path, "constraints/input_permissive_exception.yaml")
SQUARE_PLANAR_ARGS = os.path.join(test_path, "constraints/input_square_planar.yaml")
TETRAHEDRAL_ARGS = os.path.join(test_path, "constraints/input_tetrahedral.yaml")
K_ARGS = os.path.join(test_path, "constraints/input_k.yaml")
POLARISATION_ARGS = os.path.join(test_path, "constraints/input_square_planar_polarisation.yaml")
IGNORE_ARGS = os.path.join(test_path, "constraints/input_ignore.yaml")

IGNORE = "A:2002:MG__"

PASS_METAL_CONSTR = [
        '.+ 2.72\\d+,.+\"A:40:_OG_\", .+\"A:2002:MG__\".+',
        '.+ 1.99\\d+,.+\"Z:2001:_O5_\", .+\"A:2002:MG__\".+',
        '.+ 2.10\\d+,.+\"Z:2001:_O1_\", .+\"A:2002:MG__\".+',
        '.+ 2.23\\d+,.+\"A:17:_OG1", .+\"A:2002:MG__\".+'
]

METAL_CONSTR = [
       '.+ 2.03\d+,.+\"A:239:_OD1\".+\"A:350:MG__\".+',
       '.+ 2.04\d+,.+\"A:311:_OW_\".+\"A:350:MG__\".+',
       '.+ 2.53\d+,.+\"A:401:CL__\", .+\"A:350:MG__\".+',
       '.+ 2.05\d+,.+\"A:312:_OW_\", .+\"A:350:MG__\".+',
       '.+ 2.09\d+,.+\"A:141:_OG_\", .+\"A:350:MG__\".+',
       '.+ 2.09\d+,.+\"A:139:_OG_\", .+\"A:350:MG__\".+'
]

ALL_METAL_CONSTR = [
        '.+ 1.92\\d+,.+\"A:268:_NE2", .+\"A:511:ZN__\".+',
        '.+ 2.28\\d+,.+\"A:609:_OW_\", .+\"A:512:ZN__\".+',
        '.+ 2.44\\d+,.+\"A:435:_NE2", .+\"A:512:ZN__\".+',
        '.+ 2.74\\d+,.+\"A:766:_OW_\", .+\"A:512:ZN__\".+',
        '.+ 2.05\\d+,.+\"A:294:_OE1", .+\"A:511:ZN__\".+'
]

SQUARE_PLANAR = [
        '.+ 2.74\\d+,.+\"A:107:_OD2", .+\"A:302:MG__\".+',
        '.+ 2.58\\d+,.+\"A:301:_O2G", .+\"A:302:MG__\".+',
        '.+ 2.18\\d+,.+\"A:546:_OW_\", .+\"A:302:MG__\".+',
        '.+ 2.66\\d+,.+\"A:301:_O1B", .+\"A:302:MG__\".+'
]

TETRAHEDRAL = [
        '.+ 2.19\\d+,.+\"A:1081:_SG_\", .+\"A:1201:ZN__\".+',
        '.+ 2.27\\d+,.+\"A:1089:_SG_\", .+\"A:1201:ZN__\".+',
        '.+ 2.31\\d+,.+\"A:1092:_SG_\", .+\"A:1201:ZN__\".+',
        '.+ 2.27\\d+,.+\"A:1084:_ND1", .+\"A:1201:ZN__\".+'
]

K_CONSTR = [
        '.+ 1.99\\d+,.+\"B:709:_OW_\", .+\"B:603:_K__\".+',
        '.+ 2.00\\d+,.+\"B:701:_OW_\", .+\"B:603:_K__\".+',
        '.+ 2.58\\d+,.+\"B:177:_OD2", .+\"B:603:_K__\".+',
        '.+ 2.01\\d+,.+\"B:713:_OW_\", .+\"B:603:_K__\".+',
        '.+ 2.65\\d+,.+\"B:153:_OD2", .+\"B:603:_K__\".+'
]

POLARISATION = ["    1   1.6445   0.8750  0.200000 0.9545   0.8222   0.005000000   0.000000000"]

def test_metal_constraints(ext_args=METAL_CONSTR_ARGS):
    # checks metal constraints without any flags
    errors = []
    job, _ = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", METAL_CONSTR, errors)
    assert not errors


def test_no_metal_constraints(ext_args=NO_METAL_CONSTR_ARGS):
    # checks no_metal_constraints flag
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", METAL_CONSTR, errors)
    assert errors 


def test_permissive_constraints(passed=PASS_PERMISSIVE_METAL_CONSTR_ARGS, failed=FAIL_PERMISSIVE_METAL_CONSTR_ARGS):
    
    # should add constraints around the metal
    errors = []
    job = main.run_platform_from_yaml(passed)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", PASS_METAL_CONSTR, errors)
    assert not errors


def test_all_metal_constraints(ext_args=ALL_METAL_CONSTR_ARGS, ext_args_permissive=PERMISSIVE_EXCEPTION):

    # checks constrain_all_metals -> should add whatever atoms in range
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", ALL_METAL_CONSTR, errors)
    assert not errors

    # same system, but permissive -> should fail due to lack of geometry
    try:
        job = main.run_platform_from_yaml(ext_args_permissive)
    except ce.NoGeometryAroundMetal:
        assert ce.NoGeometryAroundMetal
        return
    assert False


def test_square_planar(ext_args=SQUARE_PLANAR_ARGS):
                                                                                                                                               
    errors = []                                                                                                                                                                         
    job = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", SQUARE_PLANAR, errors)
    assert not errors 


def test_tetrahedral(ext_args=TETRAHEDRAL_ARGS):                                                                                                                                    
                                                                                                                                                                                            
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    errors = tk.check_file_regex(job.pele_dir, "pele.conf", TETRAHEDRAL, errors)
    assert not errors


def test_ignore_external(ext_args=IGNORE_ARGS):

    metal_lines = []
    job = main.run_platform_from_yaml(ext_args)
    path = os.path.join(job.pele_dir, "pele.conf")
    
    with open(path, "r") as file:
        lines = file.readlines()

        for line in lines:
            if IGNORE in line:
                metal_lines.append(line)
    assert len(metal_lines) == 1


def test_polarisation(ext_args_true=POLARISATION_ARGS, ext_args_false=SQUARE_PLANAR_ARGS):

    # no polarisation
    job1 = main.run_platform_from_yaml(ext_args_false)
    mg_template_file_false = glob.glob(os.path.join(job1.pele_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz"))
    assert not mg_template_file_false

    # polarisation with factor 10
    errors = []
    job2 = main.run_platform_from_yaml(ext_args_true)
    errors = tk.check_file_regex(job2.pele_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz", POLARISATION, errors)
    assert not errors
