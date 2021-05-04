import os
import shutil
import subprocess
import warnings

from Bio.PDB import PDBParser

from pele_platform.Errors import custom_errors
from pele_platform.constants import constants


class Parameterizer:
    # Fallback paths, if only external templates available
    OPLS_IMPACT_TEMPLATE_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
    OFF_IMPACT_TEMPLATE_PATH = "DataLocal/Templates/OpenFF/Parsley/"
    ROTAMER_LIBRARY_PATH = "DataLocal/LigandRotamerLibs/"
    OBC_TEMPLATE_PATH = "DataLocal/OBC/ligandParams.txt"

    # Available methods of charge parametrization
    charge_parametrization_methods = ["am1bcc", "gasteiger", "opls2005"]

    def __init__(self, pdb_file, forcefield="OPLS2005",
                 charge_parametrization_method="am1bcc",
                 gridres=10, solvent=None, external_templates=None,
                 external_rotamers=None, as_datalocal=False,
                 pele_dir=None, exclude_terminal_rotamers=True,
                 ligand_core_constraints=None, ligand_resname=None, ligands_to_skip=None):
        """
        Initializes Parametrization to generate template and rotamer files.

        Parameters
        ----------
        pdb_file : str
            Path to the PDB file from which all HET groups will be extracted
            and parametrized
        forcefield : str
            User-defined string to select forcefield for parametrization.
            Default is "OPLS2005"
        charge_parametrization_method : str
            User-defined string to select charge parametrization method.
            Default is "am1bcc"
        gridres : int
            Resolution of the rotamers when sampling. Default is 10 degrees
        solvent : str
            Simulation solvent. Default is None and it will use
            default values: "OBC" if using any OpenFF force field and
            "VDGBNP" otherwise
        external_rotamers : list[str]
            List of paths to external rotamer files. Default is None
        external_templates : list[str]
            List of paths to external template files. Default is None
        as_datalocal : bool
            Save output files to DataLocal folder. Default is False
        pele_dir : str
            Path to PELE directory, e.g. LIG_Pele. Default is None
        exclude_terminal_rotamers : bool
            Toggle to exclude terminal rotamers. Default is True
        ligand_core_constraints : list[str]
            List of PDB atom names to be constrained as core. Default is None
        ligand_resname : str
            Residue name of the ligand. Default is None
        ligands_to_skip : List[str]
            List of residue names to skip.
        """
        self.pdb = pdb_file
        self.check_system_protonation()
        self.forcefield = self._retrieve_forcefield(forcefield)
        self.charge_parametrization_method = \
            self._check_charge_parametrization_method(
                charge_parametrization_method,
                forcefield)
        self.gridres = gridres
        self.solvent = self._retrieve_solvent_model(solvent, forcefield)

        if external_templates is not None:
            self.external_templates = external_templates
        else:
            self.external_templates = list()

        if external_rotamers is not None:
            self.external_rotamers = external_rotamers
        else:
            self.external_rotamers = list()

        self.as_datalocal = as_datalocal
        self.ligand_core_constraints = \
            self._fix_atom_names(ligand_resname, ligand_core_constraints,
                                 self.pdb)
        self.hetero_molecules = self.extract_ligands(
            pdb_file=self.pdb,
            gridres=self.gridres,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            ligand_resname=ligand_resname,
            ligand_core_constraints=self.ligand_core_constraints)

        self.pele_dir = pele_dir
        self.exclude_terminal_rotamers = exclude_terminal_rotamers
        self.ligands_to_skip = [ligand.upper() for ligand in ligands_to_skip] if ligands_to_skip else []

    @classmethod
    def from_parameters(cls, parameters):
        """
        Initializes Parametrization from simulation parameters (e.g. passed
        from Adaptive.simulation).

        Parameters
        ----------
        parameters : ParametersBuilder object
            Simulation parameters passed from Adaptive.simulation

        Returns
        -------
        obj : Parametrization object
            Parametrization object initialized from simulation parameters
        """
        obj = Parameterizer(
            pdb_file=parameters.system,
            forcefield=parameters.forcefield,
            charge_parametrization_method=parameters.charge_parametrization_method,
            gridres=parameters.gridres,
            solvent=parameters.solvent,
            external_templates=parameters.external_template,
            external_rotamers=parameters.external_rotamers,
            as_datalocal=True,
            pele_dir=parameters.pele_dir,
            exclude_terminal_rotamers=parameters.exclude_terminal_rotamers,
            ligand_core_constraints=parameters.core,
            ligand_resname=parameters.residue,
            ligands_to_skip=parameters.skip_ligand_prep,
        )

        return obj

    @staticmethod
    def extract_ligands(pdb_file, gridres, exclude_terminal_rotamers=True,
                        ligand_core_constraints=None, ligand_resname=None):
        """
        Extracts all hetero molecules in PDB and returns them as
        peleffy.topology.Molecule objects.

        Parameters
        ----------
        pdb_file : str
            Path to PDB file
        gridres : int
            Resolution of the rotamers when sampling
        exclude_terminal_rotamers : bool
            Toggle to exclude terminal rotamers. Default is True
        ligand_core_constraints : list[str]
            List of PDB atom names to be constrained as core. Default is None
        ligand_resname : str
            Residue name of the ligand. Default is None

        Returns
        -------
        unique_molecules : list[peleffy.topology.Molecule]
            List of hetero molecules extracted from the PDB file,
            without any duplicates, water molecules or single atom
            anions (e.g. Cl-, F-, etc.)
        """
        from peleffy.utils.input import PDBFile

        if not ligand_core_constraints:
            ligand_resname = None

        reader = PDBFile(pdb_file)
        molecules = reader.get_hetero_molecules(
            rotamer_resolution=gridres,
            allow_undefined_stereo=True,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            ligand_core_constraints=ligand_core_constraints,
            ligand_resname=ligand_resname)

        # Filter out water molecules and single ions.
        to_remove = []
        for molecule in molecules:
            if molecule.tag == "HOH":
                to_remove.append(molecule)
            if (molecule.tag.upper().strip() in constants.ions9
                    and len(molecule.get_pdb_atom_names()) == 1):
                to_remove.append(molecule)

        molecules = [molecule for molecule in molecules
                     if molecule not in to_remove]

        # Remove duplicates according to their molecule tag
        unique_molecules = []

        for molecule in molecules:
            if molecule.tag not in [unique.tag for unique in unique_molecules]:
                unique_molecules.append(molecule)

        return unique_molecules

    def _copy_external_parameters(self, rotamer_path, template_paths):
        """
        Copy external ligand templates and rotamers specified by the user
        to the directory that matches with the force field we work with.

        Parameters
        ----------
        rotamer_path : str
            Path to rotamers directory in pele_dir.
        template_paths : List[str]
            Paths to template directories in pele_dir (OPLS2005 and OpenFF)

        Raises
        ------
        TemplateFileNotFound if template file is not found in the path
            supplied by the user
        RotamersFileNotFound if rotamer library file is not found in the
            path supplied by the user
        """
        print("Copying external template files:",
              ", ".join(self.external_templates))

        if self.external_templates:
            for file in self.external_templates:
                try:
                    if "opls2005" in self.forcefield.type.lower():
                        shutil.copy(file, template_paths[0])
                    elif "openFF" in self.forcefield.type.lower():
                        shutil.copy(file, template_paths[-1])
                    # Copy into both OPLS2005 and OpenFF template directories,
                    # since we don't know which force field we are supposed
                    # to use
                    else:
                        shutil.copy(file, template_paths[0])
                        shutil.copy(file, template_paths[-1])
                except IOError:
                    raise custom_errors.TemplateFileNotFound(
                        f"Could not locate {file} file. "
                        f"Please double-check the path.")

        if self.external_rotamers:
            for file in self.external_rotamers:
                try:
                    shutil.copy(file, rotamer_path)
                    print("Copied external rotamer files:",
                          ", ".join(self.external_rotamers))
                except IOError:
                    raise custom_errors.RotamersFileNotFound(
                        f"Could not locate {file} file. "
                        f"Please double-check the path.")

    @staticmethod
    def _check_solvent(solvent, forcefield):
        """
        Checks if solvent is compatible with the forcefield. OpenFF
        forcefield supports OBC solvent only.

        Parameters
        -----------
        solvent : str
            Solvent selected by the user
        forcefield : str
            Forcefield selected by the user

        Returns
        --------
        solvent : str
            Solvent string, if it is compatible

        Raises
        ------
        ValueError if the chosen solvent is not compatible with the
            current force field
        """
        if forcefield.lower() != "opls2005" and solvent.lower() == "vdgbnp":
            raise ValueError("OpenFF supports OBC solvent only. Change"
                             "forcefield to 'OPLS2005' or solvent to 'OBC2'.")
        else:
            return solvent.lower()

    def _retrieve_forcefield(self, forcefield_name):
        """
        Maps forcefield YAML argument with peleffy classes.

        Parameters
        ----------
        forcefield_name : str
            The force field name defined by the user in args.forcefield

        Returns
        -------
        forcefield_obj : a peleffy.forcefield object
            The corresponding forcefield object from peleffy selected by
            the user
        """
        from peleffy.forcefield import ForceFieldSelector
        from peleffy.forcefield import OPLS2005ForceField

        # If OpenFF extension is missing, add it
        if ('openff' in forcefield_name.lower() and
                not forcefield_name.lower().endswith('offxml')):
            forcefield_name += '.offxml'

        # Select force field by name
        selector = ForceFieldSelector()
        try:
            forcefield_obj = selector.get_by_name(forcefield_name)
        except ValueError:
            print(f"Warning, invalid force field supplied, using the "
                  f"default one: \'OPLS2005\'")
            forcefield_obj = OPLS2005ForceField()

        return forcefield_obj

    def _check_charge_parametrization_method(self, method, forcefield):
        """
        Checks if charge parametrization method selected by the user
        is supported.

        Parameters
        ----------
        method : str
            Method of charge parametrization selected by the user
        forcefield : str
            Forcefield selected by the user

        Returns
        -------
        method : str
            Method of charge parametrization compatible with the forcefield

        Raises
        ------
        ValueError if the selected charge parameterization method is not
            supported
        """
        # If no method is supplied, peleffy will use its default (which
        # depends on the force field employed)
        if method is None:
            return None

        if method.lower() not in self.charge_parametrization_methods:
            raise ValueError(f"Invalid charge parametrization method, "
                             f"choose one of: "
                             f"{self.charge_parametrization_methods}.")

        return method.lower()

    def _check_external_files(self):
        """
        Checks if any of the hetero molecules extracted from PDB has a
        user-defined rotamers or template file, so they can be skipped.

        Returns
        --------
        rotamers_to_skip : list[str]
            List of hetero molecules for which the rotamers have been
            supplied in an external file.
        templates_to_skip : list[str]
            List of hetero molecules for which the templates have been
            supplied in an external file.
        """
        ligands = [ligand.tag.strip() for ligand in self.hetero_molecules]

        if self.external_rotamers:
            # Get residue names from rotamer files.
            external_rotamer_residues = \
                [os.path.basename(file).split(".")[0]
                 for file in self.external_rotamers]
            rotamers_to_skip = \
                [residue for residue in external_rotamer_residues
                 if residue in ligands]
        else:
            rotamers_to_skip = list()

        if self.external_templates:
            all_templates = self.external_templates + constants.in_pele_data
        else:
            all_templates = constants.in_pele_data

        external_template_residues = \
            [os.path.basename(file).rstrip("z").upper()
             for file in all_templates]
        templates_to_skip = \
            [residue for residue in external_template_residues
             if residue in ligands]

        external_template_residues.extend(self.ligands_to_skip)
        templates_to_skip.extend(self.ligands_to_skip)

        # Remove any duplicates
        return list(set(rotamers_to_skip)), list(set(templates_to_skip))

    def generate_ligand_parameters(self):
        """
        Generates forcefield templates and rotamer files for ligands,
        then copies the one provided by the user.

        Raises
        ------
        LigandPreparationError if any error is obtained throughout
            the ligand preparation process
        """
        from peleffy.topology import RotamerLibrary, Topology
        from peleffy.utils import OutputPathHandler
        from peleffy.template import Impact
        from peleffy.solvent import OBC2, OPLSOBC

        rotamer_library_path, impact_template_paths = None, None

        rotamers_to_skip, templates_to_skip = self._check_external_files()
        topologies = list()

        for molecule in self.hetero_molecules:
            output_handler = OutputPathHandler(molecule,
                                               self.forcefield,
                                               as_datalocal=True,
                                               output_path=self.pele_dir)

            rotamer_library_path = output_handler.get_rotamer_library_path()
            impact_template_path = output_handler.get_impact_template_path()
            solvent_template_path = output_handler.get_solvent_template_path()

            if molecule.tag.strip() not in rotamers_to_skip:
                rotamer_library = RotamerLibrary(molecule)
                rotamer_library.to_file(rotamer_library_path)

            # This boolean indicates whether we needed to reparameterize
            # this ligand with OPLS2005 or not
            opls_reparameterization = False

            if molecule.tag.strip() not in templates_to_skip:
                # Try to parametrize with OPLS if OpenFF fails
                try:
                    parameters = self.forcefield.parameterize(
                        molecule,
                        self.charge_parametrization_method)
                except (subprocess.CalledProcessError, TypeError) as e1:
                    warnings.warn(f"Could not parameterize residue "
                                  f"{molecule.tag.strip()} with the selected "
                                  f"forcefield. The following error was "
                                  f"obtained: {e1}")

                    default = "OPLS2005"
                    fallback_forcefield = self._retrieve_forcefield(default)

                    try:
                        parameters = fallback_forcefield.parameterize(molecule)
                        warnings.warn(f"Parametrized with {default} "
                                      f"instead.")
                    except subprocess.CalledProcessError as e2:
                        raise custom_errors.LigandPreparationError(
                            f"Could not parametrize {molecule.tag.strip()}. "
                            f"The error was {e2}.")

                    opls_reparameterization = True

                topology = Topology(molecule, parameters)
                topologies.append(topology)

                # In case we reparameterized the templates with OPLS,
                # we need to convert atom types to OFFT (unique atom
                # type for OpenFF). Otherwise, PELE will complain about it
                if opls_reparameterization:
                    for atom in topology.atoms:
                        atom.set_OPLS_type('OFFT')

                impact = Impact(topology)
                impact.to_file(impact_template_path)

                impact_template_paths = [impact_template_path]

                print(f"Parametrized {molecule.tag.strip()}.")

        if self.solvent is not None:
            if self.solvent.upper() == 'OBC':
                if 'openff' in self.forcefield.type.lower():
                    solvent_parameters = OBC2(topologies)
                else:
                    solvent_parameters = OPLSOBC(topologies)

                solvent_parameters.to_file(os.path.join(self.pele_dir,
                                                        self.OBC_TEMPLATE_PATH))
        else:
            if 'openff' in self.forcefield.type.lower():
                solvent_parameters = OBC2(topologies)
                solvent_parameters.to_file(os.path.join(self.pele_dir,
                                                        self.OBC_TEMPLATE_PATH))

        if not rotamer_library_path:
            rotamer_library_path = os.path.join(self.pele_dir,
                                                self.ROTAMER_LIBRARY_PATH)

        if not impact_template_paths:
            impact_template_paths = [
                os.path.join(self.pele_dir, self.OPLS_IMPACT_TEMPLATE_PATH),
                os.path.join(self.pele_dir, self.OFF_IMPACT_TEMPLATE_PATH)]

        self._copy_external_parameters(os.path.dirname(rotamer_library_path),
                                       [os.path.dirname(path)
                                        for path in impact_template_paths])

    def _retrieve_solvent_model(self, solvent_name, forcefield):
        """
        Checks solvent compatibility with the forcefield and returns the
        solvent class from peleffy.

        Parameters
        -----------
        solvent_name : str
            Name of solvent defined by the user in YAML
        forcefield : str
            Forcefield defined by the user in YAML

        Returns
        --------
        solvent_class : a peleffy.solvent object
            The solvent object from peleffy that is selected
        """
        from peleffy.solvent import OBC2, OPLSOBC

        solvent_class = None

        if solvent_name is None:
            if 'openff' in forcefield.lower():
                solvent_class = OBC2
        else:
            checked_solvent = self._check_solvent(solvent_name, forcefield)

        if forcefield.lower() == 'opls2005' and solvent == 'obc':
            solvent_class = OPLSOBC

        return solvent_class

    @staticmethod
    def _fix_atom_names(ligand_resname, ligand_core_constraints, pdb):
        """
        Adds spaces around PDB atom names to ensure peleffy can correctly
        identify core atoms.

        Parameters
        -----------
        ligand_resname : str
            Residue name of the ligand
        ligand_core_constraints : list[str]
            List of PDB atom names to constrain as core defined by the
            user, e.g. ["N1", "O1"]
        pdb : str
            Path to PDB file (self.pdb)

        Returns
        --------
        ligand_core_constraints : list[str]
            List of fixed PDB atom names of the ligand

        Raises
        ------
        ValueError if any constraint points to an atom missing in the PDB
        """
        if ligand_core_constraints:
            parser = PDBParser()
            structure = parser.get_structure("system", pdb)
            user_constraints = [atom.strip()
                                for atom in ligand_core_constraints]
            fixed_atoms = []

            for residue in structure.get_residues():
                if residue.get_resname() == ligand_resname:
                    for atom in residue.get_atoms():
                        if atom.name in user_constraints:
                            fixed_atoms.append(atom.fullname)

            # If the number of fixed atoms doesn't match the input,
            # raise an Error with a list of missing atoms.
            if len(user_constraints) != len(fixed_atoms):
                not_found = [atom for atom in user_constraints
                             if atom not in [pdb_atom.strip()
                                             for pdb_atom in fixed_atoms]]
                raise ValueError(f"Atom(s) {not_found} were not "
                                 f"found in {pdb}.")

            return fixed_atoms

    def check_system_protonation(self):
        """
        It checks for hydrogen atoms in the input PDB to ensure that
        it has been protonated and complaints otherwise.

        Raises
        ------
        ProtonationError if no hydrogen atoms are found in the input PDB
        """
        with open(self.pdb) as f:
            lines = f.readlines()
            atom_lines = [line for line in lines
                          if line.startswith("ATOM")
                          or line.startswith("HETATM")]

            for line in atom_lines:
                if line[12:16].strip().startswith("H"):
                    return

        raise custom_errors.ProtonationError("We did not find any hydrogens "
                                             "in your system - looks like "
                                             "you forgot to "
                                             "protonate it.")
