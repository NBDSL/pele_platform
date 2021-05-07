import os
import tempfile
import shutil
import sys

import pele_platform.Utilities.Helpers.simulation as ad
import pele_platform.Frag.helpers as hp
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce
import pele_platform.Frag.libraries as lb
import pele_platform.Frag.analysis as ana
import frag_pele.main as frag
import pele_platform.constants.constants as cs


class FragRunner(object):

    def __init__(self, args):
        from pele_platform.Utilities.Parameters import ParametersBuilder

        builder = ParametersBuilder()
        parameters = builder.build_frag_variables(args)
        self._parameters = parameters

    @property
    def parameters(self):
        return self._parameters

    def run_simulation(self):
        self._set_test_variables()
        self._prepare_control_file()
        self._launch()
        self._analysis()

        return self.parameters

    def _launch(self):
        params = self.parameters

        with tempfile.TemporaryDirectory() as tmpdirname:
            if params.ligands:  # Full ligands as sdf
                fragment_files = self._prepare_input_file(logger=params.logger)
            elif params.frag_library:
                params.input = lb.main(params.core,
                                       params.frag_core_atom,
                                       params.frag_library,
                                       params.logger,
                                       params.fragment_atom,
                                       params.frag_restart,
                                       tmpdirname)
            else:
                fragment_files = None
                assert params.input is not None, "You need to provide input.conf file or SD file with fully grown " \
                                                 "ligands. "

            if not params.only_analysis:
                self._run()

            if params.cleanup and fragment_files:
                self._clean_up(fragment_files)

    def _prepare_control_file(self):
        # Create tmp folder with frag control_file
        params = self.parameters

        tmp_dir = tempfile.mkdtemp()
        tmp_control_file = os.path.join(tmp_dir,
                                        os.path.basename(params.control_file))
        shutil.copy(params.control_file, tmp_control_file)
        adaptive = ad.SimulationBuilder("", tmp_control_file, params)
        # Fill to time because we have flags inside flags
        adaptive.fill_pele_template(params, params.water_object)
        adaptive.fill_pele_template(params, params.water_object)
        self.parameters.control_file = tmp_control_file

        return self.parameters.control_file

    def _set_test_variables(self):
        if self.parameters.test:
            self.parameters.gr_steps = 1
            self.parameters.frag_steps = 2
            self.parameters.frag_eq_steps = 1
            self.parameters.temperature = 100000
            self.parameters.anm_freq = 0
            self.parameters.minimizatoon = 0
            self.parameters.sidechain_freq = 0
            self.parameters.water_freq = 0
            self.parameters.cpus = 4
        else:
            pass

    def _run(self):
        params = self.parameters
        params.spython = cs.SCHRODINGER
        if params.frag_library:
            params.working_dir = []
            pdb_basename = params.core.split(".pdb")[0] + "_processed" # Get the name of the pdb without extension
            if "/" in pdb_basename:
                pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name
            current_path = os.path.abspath(".")
            with open(params.input, "r") as input_file:
                for line in input_file.readlines():
                    ID = line.split('/')[-1].split(".pdb")
                    ID = "".join(ID[:]).replace(" ", "")
                    params.working_dir.append(os.path.join(current_path, "{}_{}".format(pdb_basename, ID)).strip('\n'))
        else:
            pdb_basename = params.core.split(".pdb")[0] + "_processed"  # Get the name of the pdb without extension
            if "/" in pdb_basename:
                pdb_basename = pdb_basename.split("/")[-1] # And if it is a path, get only the name
            current_path = os.path.abspath(".")
            ID = open(params.input,'r').readlines()[0].split('.pdb')
            ID = "".join(ID[:]).replace(" ","")
            params.working_dir = os.path.join(current_path, "{}_{}".format(pdb_basename, ID)).strip('\n')
        if params.frag_run:
            try:
                frag.main(params.core_process, params.input, params.gr_steps,
                          params.criteria, params.plop_path, params.spython,
                          params.pele_exec, params.control_file,
                          params.license, params.output_folder,
                          params.report_name, "trajectory",
                          params.cluster_folder, params.cpus, params.distcont,
                          params.threshold, params.epsilon, params.condition,
                          params.metricweights, params.nclusters,
                          params.frag_eq_steps, params.frag_restart,
                          params.min_overlap, params.max_overlap,
                          params.chain_core, params.frag_chain,
                          params.frag_steps, params.temperature, params.seed,
                          params.gridres, params.banned, params.limit,
                          params.mae, params.rename, params.threshold_clash,
                          params.steering, params.translation_high,
                          params.rotation_high, params.translation_low,
                          params.rotation_low, params.explorative,
                          params.frag_radius, params.sampling_control,
                          params.pele_data, params.pele_documents,
                          params.only_prepare, params.only_grow,
                          params.no_check, params.debug, srun=params.usesrun)
            except Exception as e:
                print("Skipped - FragPELE will not run.")
                print(e)

    def _prepare_input_file(self, logger=None):
        from rdkit import Chem
        import rdkit.Chem.rdmolops as rd

        self.parameters.input = "input.conf"

        params = self.parameters

        # Check input file
        limit_atoms_ligand = 100
        ch.check_limit_number_atoms(params.ligands, limit_atoms_ligand)

        # Get core of the ligand
        mol = Chem.MolFromPDBFile(params.core)
        try:
            ligand_core = rd.SplitMolByPDBResidues(mol)[params.residue]
        except KeyError:
            raise ce.MissResidueFlag("Missing residue flag to specify " +
                                     "the ligand core residue name. " +
                                     "i.e resname: 'LIG'")

        # Get sdf full grown ligands
        ligands_grown = Chem.SDMolSupplier(params.ligands, removeHs=False)
        fragment_files = []

        # For each full grown ligand create neutral fragment
        with open(params.input, "w") as fout:
            pass
        for ligand in ligands_grown:
            Chem.AssignAtomChiralTagsFromStructure(ligand)
            try:
                line, fragment = \
                    self._create_fragment_from_ligand(ligand,
                                                      ligand_core)
            except Exception as e:
                try:
                    line, fragment = \
                        self._create_fragment_from_ligand(ligand,
                                                          ligand_core,
                                                          substructure=False)
                except Exception as e:
                    try:
                        # Try to fix symmetry
                        line, fragment = \
                            self._create_fragment_from_ligand(ligand,
                                                              ligand_core,
                                                              symmetry=True)
                    except Exception as e:
                        # Try with second substructure search
                        line, fragment = \
                            self._create_fragment_from_ligand(
                                ligand, ligand_core,
                                result=1, substructure=False)

            logger.info(f"Ligand {fragment.file} preprocessed")
            fragment_files.append(fragment.file)
            with open(params.input, "a") as fout:
                fout.write(line + "\n")

        return fragment_files

    def _create_fragment_from_ligand(self, ligand, ligand_core, result=0, substructure=True, symmetry=False):
        import rdkit.Chem.AllChem as rp
        from rdkit import Chem

        params = self.parameters

        fragment, old_atoms, hydrogen_core, atom_core, atom_frag, \
        mapping, correct = hp._build_fragment_from_complex(
            params.core, params.residue, ligand, ligand_core,
            result, substructure, symmetry, params.frag_core_atom)

        # temporary override to fix segmentation faults
        filename = "temp.pdb"
        Chem.MolToPDBFile(fragment, filename)
        fragment = Chem.MolFromPDBFile(filename, removeHs=False)
        os.remove(filename)

        rp.EmbedMolecule(fragment)
        fragment = hp._retrieve_fragment(
            fragment, old_atoms, atom_core, hydrogen_core, atom_frag, mapping)
        line = fragment.get_inputfile_line()
        fragment.sanitize_file()

        if not correct:
            print("Ligand incorrect")
        return line, fragment

    def _analysis(self):
        def run_analysis(sim_directory, simulation_output):

            # Only proceed if sampling_result folder exists


            analysis_folder = os.path.join(sim_directory, "results")

            analysis = Analysis(resname='GRW',
                                chain=self.parameters.chain,
                                simulation_output=simulation_output,
                                be_column=self.parameters.be_column,
                                limit_column=self.parameters.limit_column,
                                traj=self.parameters.traj_name,
                                report=self.parameters.report_name,
                                skip_initial_structures=not self.parameters.test,
                                kde=self.parameters.kde,
                                kde_structs=self.parameters.kde_structs,
                                topology=self.parameters.topology,
                                cpus=self.parameters.cpus)
            analysis.generate(
                analysis_folder,
                clustering_type=self.parameters.clustering_method.lower(),
                bandwidth=self.parameters.bandwidth,
                analysis_nclust=self.parameters.analysis_nclust,
                max_top_clusters=self.parameters.max_top_clusters,
                top_clusters_criterion=self.parameters.top_clusters_criterion,
                min_population=self.parameters.min_population,
                representatives_criterion=self.parameters.cluster_representatives_criterion)

        self.parameters.analysis_to_point = \
            self.parameters.args.analysis_to_point

        if self.parameters.analysis_to_point and self.parameters.folder:
            ana.main(path=self.parameters.folder,
                     atom_coords=self.parameters.analysis_to_point,
                     pattern=os.path.basename(self.parameters.system))

        # TODO create a list of the libraries defined in the current input.yaml
        from pele_platform.analysis import Analysis
        from glob import glob
        if self.parameters.frag_library:
            for sim_directory in self.parameters.working_dir:
                simulation_output = sim_directory + '/sampling_result'
                run_analysis(sim_directory, simulation_output)
        else:
            simulation_output = self.parameters.working_dir +  '/sampling_result'
            if not os.path.isdir(simulation_output):
                return 0
            run_analysis(self.parameters.working_dir,simulation_output)


    def _clean_up(self, fragment_files):
        for file in fragment_files:
            if os.path.isfile(file):
                os.remove(file)
