import os
import pytest

from pele_platform.Adaptive import parameterizer
from pele_platform.constants import constants
from pele_platform import main

test_path = os.path.join(constants.DIR, "Examples", "constraints")


@pytest.mark.parametrize(
    ("pdb", "ligands", "core", "resname"),
    [
        (
            os.path.join(test_path, "1zop.pdb"),
            ["MG", "LIG"],
            [" O1 ", " N1 ", " C1 ", " C2 "],
            "LIG",
        ),
        (
            os.path.join(test_path, "4zu9_prep.pdb"),
            ["MET", "MG", "SO4", "LIG"],
            None,
            None,
        ),
        (
            os.path.join(test_path, "4zu9_prep.pdb"),
            ["MET", "MG", "SO4", "LIG"],
            None,
            "LIG",
        ),
    ],
)
def test_ligand_extraction(pdb, ligands, core, resname):
    """
    Checks hetatm ligand extraction from PDB file in Parametrization class.

    Parameters
    -----------
    pdb : str
        Path to PDB file.
    ligands : List[str]
        List of expected residue names to be extracted from the PDB.
    """
    extracted_ligands = parameterizer.Parameterizer.extract_ligands(
        pdb, gridres=10, ligand_core_constraints=core, ligand_resname=resname
    )
    assert len(extracted_ligands) == len(ligands)


def test_generate_ligand_parameters(parametrize):
    """
    Runs the whole parametrization flow and checks if all expected files were created.

    Parameters
    -----------
    parametrize : Parametrization
        Object created in fixture.
    """
    expected_templates = ["ligz"]  # MG template included in pele Data
    expected_rotamers = ["MG.rot.assign", "LIG.rot.assign"]
    parametrize.generate_ligand_parameters()

    # Check if all expected files were created and clean them up
    for file in expected_templates + expected_rotamers:
        assert os.path.exists(file)
        os.remove(file)


@pytest.mark.parametrize(
    ("solvent", "forcefield", "error"),
    [
        ("OBC2", "OPLS2005", False),
        ("vdgbnp", "OPLS2005", False),
        ("OBC2", "openff-1.2.0", False),
        ("vdgbnp", "openff-1.3.0", True),
    ],
)
def test_check_solvent(solvent, forcefield, error):
    """
    Checks if certain combinations of forcefield and solvent raise error when necessary.

    Parameters
    -----------
    solvent : str
        User-defined solvent, one of parametrization.forcefields.
    forcefield : str
        User-defined solvent, one of parametrization.solvents
    error : bool
        Set to True if particular combination is supposed to raise en error.
    """
    if error:
        with pytest.raises(ValueError):
            parameterizer.Parameterizer._check_solvent(
                solvent=solvent, forcefield=forcefield
            )
    else:
        parameterizer.Parameterizer._check_solvent(
            solvent=solvent, forcefield=forcefield
        )


@pytest.mark.parametrize(
    ("method", "forcefield", "error", "expected_output"),
    [
        ("BLABLABLA", "OPLS2005", True, None),
    ],
)
def test_check_charge_parametrization_method(
    parametrize, method, forcefield, error, expected_output
):
    """
    Tests the checker for charge parametrization method.

    Parameters
    ------------
    parametrize : parametrize.Parametrization
        Object created in the fixture.
    method : str
        User-defined charge parametrization method.
    forcefield : str
        User-defined forcefield
    error : bool
        True if we're expecting the combination of arguments to raise an error.
    expected_output : str
        Expected compatible method to be returned.
    """
    if error:
        with pytest.raises(ValueError):
            parametrize._check_charge_parametrization_method(method, forcefield)
    else:
        output = parametrize._check_charge_parametrization_method(method, forcefield)
        assert output == expected_output


@pytest.mark.parametrize(
    ("user_input", "expected_output", "resname", "pdb"),
    [
        (
            ["O1", "N1", "C1", "C2"],
            [" O1 ", " N1 ", " C1 ", " C2 "],
            "LIG",
            os.path.join(test_path, "1zop.pdb"),
        ),
        (
            ["C14", "C15", "O4", "O5"],
            [" C14", " C15", " O4 ", " O5 "],
            "LIG",
            os.path.join(test_path, "4l3a_prep.pdb"),
        ),
    ],
)
def test_fix_atom_names(user_input, expected_output, resname, pdb):
    """
    Tests the Parametrization method for fixing core atom names provided by the user.

    Parameters
    ----------
    user_input : List[str]
        List of atom names to constrain when parametrizing the ligand.
    expected_output : List[str]
        User input with added spaces to match peleffy specifications.
    resname : str
        Ligand residue name.
    pdb : str
        Path to PDB file.
    """
    output = parameterizer.Parameterizer._fix_atom_names(
        ligand_resname=resname, ligand_core_constraints=user_input, pdb=pdb
    )
    assert expected_output == output


def test_production():
    """
    Runs full production from YAML with all flags.
    TODO: Uncomment forcefield and solvent in YAML whenever OBC2 starts working...
    """
    yaml_path = os.path.join(
        constants.DIR, "Examples", "preparation", "parametrization_flags.yaml"
    )
    main.run_platform_from_yaml(os.path.join(yaml_path))


@pytest.fixture
def parametrize():
    """
    Creates Parametrization class to use in other tests.

    Returns
    --------
    obj : parametrization.Parameterizer
        Parametrization object to use in tests.
    """
    obj = parameterizer.Parameterizer(
        pdb_file=os.path.join(test_path, "1zop.pdb"),
        forcefield="openff-1.3.0",
        solvent="OBC",
        as_datalocal=False,
    )
    return obj
