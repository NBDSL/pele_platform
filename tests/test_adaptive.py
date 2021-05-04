import os
import math
import glob
import shutil
import re

import pele_platform.constants.constants as cs
import pele_platform.main as main

test_path = os.path.join(cs.DIR, "Examples")


OUT_IN_ARGS = os.path.join(test_path, "out_in/input.yaml")
INDUCED_EX_ARGS = os.path.join(test_path, "induced_fit/input_exhaustive.yaml")
INDUCED_FAST_ARGS = os.path.join(test_path, "induced_fit/input_fast.yaml")
NWATER_ARGS = os.path.join(test_path, "water/input_nwaters.yaml")
GLOBAL_ARGS = os.path.join(test_path, "global/input.yaml")
INPUTS_GLOBAL_ARGS = os.path.join(test_path, "global/input_inputs.yaml")
INPUTS_AST_ARGS = os.path.join(test_path, "global/input_inputs_asterisc.yaml")
EXIT_ARGS = os.path.join(test_path, "exit/input.yaml")
EXITSOFT_ARGS = os.path.join(test_path, "exit_soft/input.yaml")
WATER_ARGS = os.path.join(test_path, "water/input_bs.yaml")
ALL_WATER_ARGS = os.path.join(test_path, "water/input_all.yaml")
WATERLIG_ARGS = os.path.join(test_path, "water/input_lig.yaml")
RESTART_ARGS = os.path.join(test_path, "restart/input.yaml")
MSM_ARGS = os.path.join(test_path, "Msm/input.yaml")
MAE_ARGS = os.path.join(test_path, "induced_fit/input_mae.yaml")
PCA_ARGS = os.path.join(test_path, "pca/input.yaml")
PCA2_ARGS = os.path.join(test_path, "pca/input_str.yaml")
FLAGS_ARGS = os.path.join(test_path, "flags/input.yaml")
RESCORING_ARGS = os.path.join(test_path, "rescoring/input.yaml")
GPCR_ARGS = os.path.join(test_path, "gpcr/input.yaml")

ADAPTIVE_VALUES = [
    "water_processed.pdb",
    "SB4",
    '"outputPath": "output_sim"',
    '"processors" : 3',
    '"peleSteps" : 1,',
    '"iterations" : 1,',
    '"runEquilibration" : true,',
    '"equilibrationLength" : 11,',
    '"seed": 3000',
    '"values" : [1, 2, 3],',
    '"conditions": [0.1, 0.2, 0.3]',
    '"epsilon": 0.3',
    '"metricColumnInReport" : 3,',
    '"metricColumnInReport" : 3,',
    '"type" : "epsilon"',
    "exitContinuous",
    '"data": "done"',
    '"executable": "done"',
    '"documents": "done"',
]

PELE_VALUES = [
    "rep",
    "traj.xtc",
    "OBC",
    '"anmFrequency" : 3,',
    '"sideChainPredictionFrequency" : 3,',
    '"minimizationFrequency" : 3,',
    '"temperature": 3000,',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:111:_CA_" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:11:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 3, "equilibriumDistance": 0.0, "constrainThisAtom": "A:13:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:353:_CA_" }',
    '{ "type": "constrainAtomToPosition", "springConstant": 5.0, "equilibriumDistance": 0.0, "constrainThisAtom": "A:5:_CA_" }',
    '"radius": 3000',
    '"fixedCenter": [30,30,30]',
    'tests/native.pdb"',
    '"atoms": { "ids":["A:6:_CG_"]}',
    '"atoms": { "ids":["A:6:_CD_"]}',
    "simulationLogPath",
    '"activateProximityDetection": false',
    '"overlapFactor": 3',
    '"steeringUpdateFrequency": 3',
    '"verboseMode": true,',
    '"displacementFactor" : 3',
    '"modesChangeFrequency" : 3,',
    '"directionGeneration" : "steered",',
    '"modesMixingOption" : "DontMixModes",',
    '"pickingCase" : "random"',
    '"numberOfModes": 3',
    '"preloadedModesIn" : "modes.nmd",',
    '"relaxationSpringConstant" : 3',
]

PCA_VALUES = [
    '"displacementFactor" : 2',
    '"modesChangeFrequency" : 50,',
    '"directionGeneration" : "oscillate",',
    '"modesMixingOption" : "doNotMixModes",',
    '"pickingCase" : "LOWEST_MODE"',
    '"numberOfModes": 1',
    '"relaxationSpringConstant" : 2',
]

ALL_WATER_VALUES = ["WaterPerturbation::parameters", '"M:1"', '"M:2"']

WATER_VALUES = [
    "WaterPerturbation::parameters",
    '"M:1"',
]

GPCR_VALUES = [
    '"radius": 19.970223159033843,',
    '"fixedCenter": [-71.78435134887695,-13.431749963760375,-42.46209926605225]',
]


def test_induced_exhaustive(ext_args=INDUCED_EX_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_induced_fast(ext_args=INDUCED_FAST_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_n_water(ext_args=NWATER_ARGS):
    job = main.run_platform_from_yaml(ext_args)
    results = glob.glob(os.path.join(job.pele_dir, "results/top_poses/*.pdb"))
    error = False
    # Result has waters
    for result in results:
        with open(result, "r") as f:
            if not "HOH" in "".join(f.readlines()):
                error = True
    # Input has no water
    with open("../pele_platform/Examples/Msm/PR_1A28_xray_-_minimized.pdb", "r") as f:
        if "HOH" in "".join(f.readlines()):
            error = True
    assert not error


def test_global(ext_args=GLOBAL_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_inputs_global(ext_args=INPUTS_GLOBAL_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_exit(ext_args=EXIT_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_softexit(ext_args=EXITSOFT_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_all_waters(ext_args=ALL_WATER_ARGS):
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    folder = job.pele_dir
    errors = check_file(folder, "pele.conf", ALL_WATER_VALUES, errors)
    assert not errors


def test_water_lig(ext_args=WATERLIG_ARGS):
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    folder = job.pele_dir
    errors = check_file(folder, "pele.conf", WATER_VALUES, errors)
    assert not errors


def test_out_in(ext_args=OUT_IN_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_restart(ext_args=RESTART_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_rescoring(ext_args=RESCORING_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_mae(ext_args=MAE_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_asterisc_inputs(ext_args=INPUTS_AST_ARGS):
    main.run_platform_from_yaml(ext_args)


def test_flags(ext_args=FLAGS_ARGS, output="NOR_solvent_OBC"):
    errors = []
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    job = main.run_platform_from_yaml(ext_args)
    folder = job.folder
    if not os.path.exists(
        os.path.join(folder, "DataLocal/LigandRotamerLibs/STR.rot.assign")
    ) or not os.path.exists(
        os.path.join(folder, "DataLocal/LigandRotamerLibs/MG.rot.assign")
    ):
        errors.append("External rotamer flag not working")
    if not os.path.exists(
        os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/strz")
    ) or not os.path.exists(
        os.path.join(folder, "DataLocal/Templates/OPLS2005/HeteroAtoms/mgz")
    ):
        errors.append("External templates flag not working")
    if os.path.exists(os.path.join(folder, job.output)):
        errors.append("Debug flag not working")
    errors = check_file(folder, "adaptive.conf", ADAPTIVE_VALUES, errors)
    errors = check_file(folder, "pele.conf", PELE_VALUES, errors)
    errors = check_file(
        folder, "DataLocal/LigandRotamerLibs/SB4.rot.assign", "60", errors
    )
    assert not errors


def test_pca(ext_args=PCA_ARGS, output="PCA_result"):
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    folder = job.folder
    errors = check_file(folder, "pele.conf", PCA_VALUES, errors)
    assert not errors


def test_str_pca(ext_args=PCA2_ARGS, output="PCA_result"):
    if os.path.exists(output):
        shutil.rmtree(output, ignore_errors=True)
    errors = []
    job = main.run_platform_from_yaml(ext_args)
    folder = job.folder
    errors = check_file(folder, "pele.conf", PCA_VALUES, errors)
    assert not errors


def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n


def check_file(
    folder, filename, values, errors, subdelimiter=None, truncate_digits_to=4
):
    filename = os.path.join(folder, filename)
    print(values)
    print("---------------------")
    print("filename:", filename)
    with open(filename, "r") as f:
        lines = f.readlines()
        print("lines:", lines)
        if subdelimiter is None:
            for value in values:
                if value not in "".join(lines):
                    errors.append(value)

        elif isinstance(subdelimiter, str):
            for value in values:
                splitted_values = value.split(subdelimiter)
            for splitted_value in splitted_values:
                if splitted_value.replace(".", "").replace("-", "").isnumeric():
                    if str(
                        truncate(float(splitted_value), truncate_digits_to)
                    ) not in "".join(lines):
                        errors.append(value)
                else:
                    if splitted_value not in "".join(lines):
                        errors.append(splitted_value)
        else:
            raise TypeError("Wrong subdelimiter type")

    return errors


def test_gpcr(args=GPCR_ARGS):
    errors = []
    job = main.run_platform_from_yaml(args)
    errors = check_file(job.pele_dir, "pele.conf", GPCR_VALUES, errors)
    input_file = os.path.join(job.inputs_dir, "complex_processed.pdb")
    if not os.path.exists(input_file):
        errors.append("skip_ppp")
    assert not errors


def check_file_regex(directory, file, patterns, errors):
    """
    Checks file for strings matching regex patterns.

    Parameters
    -----------
    directory : str
        Path to folder.
    file : str
        Name of the file to check.
    patterns : List[str]
        List of regex patterns to check.
    errors : List[str]
        List of existing errors (or empty list).
    """
    path = os.path.join(directory, file)

    with open(path) as f:
        content = f.read()

    for pattern in patterns:
        if not re.search(pattern, content):
            errors.append(pattern)

    return errors
