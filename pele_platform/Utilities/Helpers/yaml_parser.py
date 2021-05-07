from dataclasses import dataclass
from pele_platform.Errors.custom_errors import LigandNameNotSupported
from difflib import SequenceMatcher
import os
import yaml
import warnings


def _yaml_error_wrapper(error):
    """
    Wraps YAML errors into a more human-friendly format, making customs suggestions about potential issues.
    This should be expanded in the future when more issues get reported by the users.
    """
    custom_errors = {
        "expected '<document start>', but found '<block mapping start>'": "Please ensure every key in input.yaml is "
        "followed by a colon and "
        "a space. There seem to be some issues on line {}, character {}.",
        "found character '\\t' that cannot start any token": "Please remove any trailing tabs from input.yaml, there "
        "seem to be one on line {}, character {}.",
    }

    custom = custom_errors.get(str(error.problem).strip(), None)

    if custom:
        line = error.problem_mark.line + 1
        character = int(error.problem_mark.column) + 1
        raise yaml.YAMLError(custom.format(line, character))
    else:
        raise error


@dataclass
class YamlParser(object):
    yamlfile: str
    valid_flags: dict

    def read(self) -> None:
        self.data = self._parse_yaml()
        self._check()
        self._check_residue()
        self._parse()
        self._get_value_from_env()

    def _parse_yaml(self) -> dict:
        # Retrieve raw info from yaml
        with open(self.yamlfile, "r") as stream:
            try:
                data = yaml.safe_load(stream)
            except Exception as error:
                _yaml_error_wrapper(error)
        return data

    def _get_value_from_env(self):
        self.usesrun = bool(os.environ.get("SRUN", self.usesrun))


    def _check(self) -> None:
        # Check if valids in yaml file are valids
        for key in self.data.keys():
            if key not in self.valid_flags.values():
                raise KeyError(self._recommend(key))

    def _check_residue(self) -> None:
        if "resname" in self.data.keys():
            if self.data["resname"] == "UNK":
                raise LigandNameNotSupported(
                    "'UNK' ligand name is not supported, please rename it, e.g. 'LIG'."
                )

    def _recommend(self, key):
        most_similar_flag = None
        for valid_key in self.valid_flags.values():
            flag = Most_Similar_Flag(valid_key)
            flag.calculate_distance(key)
            if not most_similar_flag:
                most_similar_flag = flag
            else:
                if flag.distance > most_similar_flag.distance:
                    most_similar_flag = flag
        exception_raised = (
            f"Incorrect flag {key}. Did you mean {most_similar_flag.name}?"
        )
        return exception_raised

    def _parse(self) -> None:
        # Parse fields in yaml file and set defaults
        valid_flags = self.valid_flags
        data = self.data
        self.system = data.get(valid_flags["system"], "")
        self.system = os.path.abspath(self.system) if self.system else ""
        self.residue = data.get(valid_flags["residue"], None)
        self.chain = data.get(valid_flags["chain"], None)
        self.hbond = data.get(valid_flags["hbond"], [None, None])
        self.test = data.get(valid_flags["test"], None)
        self.pele = data.get(valid_flags["pele"], None)
        self.forcefield = data.get(valid_flags["forcefield"], "OPLS2005")
        self.verbose = data.get(valid_flags["verbose"], None)
        self.anm_freq = data.get(valid_flags["anm_freq"], None)
        self.sidechain_freq = data.get(valid_flags["sidechain_freq"], None)
        self.min_freq = data.get(valid_flags["min_freq"], None)
        self.water_freq = data.get(valid_flags["water_freq"], None)
        self.temperature = self.temp = data.get(valid_flags["temperature"], None)
        self.sidechain_resolution = data.get(valid_flags["sidechain_resolution"], None)
        self.steric_trials = data.get(valid_flags["steric_trials"], None)
        self.overlap_factor = data.get(valid_flags["overlap_factor"], None)
        self.steering = data.get(valid_flags["steering"], None)
        self.solvent = data.get(valid_flags["solvent"], None)
        self.usesrun = data.get(valid_flags["usesrun"], None)
        self.spawning = data.get(valid_flags["spawning"], None)
        self.iterations = data.get(valid_flags["iterations"], None)
        self.pele_steps = self.steps = data.get(valid_flags["pele_steps"], None)
        self.cpus = data.get(valid_flags["cpus"], None)
        self.density = data.get(valid_flags["density"], None)
        self.cluster_values = data.get(valid_flags["cluster_values"], None)
        self.cluster_conditions = data.get(valid_flags["cluster_conditions"], None)
        self.simulation_type = data.get(valid_flags["simulation_type"], None)
        self.equilibration = data.get(valid_flags["equilibration"], None)
        self.clust_type = data.get(valid_flags["clust_type"], None)
        self.eq_steps = data.get(valid_flags["eq_steps"], None)
        self.adaptive_restart = data.get(valid_flags["adaptive_restart"], None)
        self.input = data.get(valid_flags["input"], None)
        self.report_name = data.get(valid_flags["report_name"], None)
        self.traj_name = data.get(valid_flags["traj_name"], None)
        self.adaptive = data.get(valid_flags["adaptive"], None)
        self.epsilon = data.get(valid_flags["epsilon"], None)
        self.out_in = data.get(valid_flags["out_in"], None)
        self.bias_column = data.get(valid_flags["bias_column"], None)
        self.gridres = data.get(valid_flags["gridres"], 10)
        self.core = data.get(valid_flags["core"], -1)
        self.mtor = data.get(valid_flags["mtor"], 4)
        self.n = data.get(valid_flags["n"], 10000)
        self.template = data.get(valid_flags["template"], None)
        self.ext_temp = self.template
        self.rotamers = data.get(valid_flags["rotamers"], None)
        self.ext_rotamers = self.rotamers
        self.mae_lig = data.get(valid_flags["mae_lig"], None)
        self.mae_lig = os.path.abspath(self.mae_lig) if self.mae_lig else None
        self.skip_prep = self.no_ppp = data.get(valid_flags["skip_prep"], None)
        self.gaps_ter = data.get(valid_flags["gaps_ter"], None)
        self.charge_ter = data.get(valid_flags["charge_ter"], None)
        self.mpi_params = data.get(valid_flags["mpi_params"], None)
        self.nonstandard = data.get(valid_flags["nonstandard"], None)
        self.prepwizard = data.get(valid_flags["prepwizard"], None)
        self.box_center = data.get(valid_flags["box_center"], None)
        self.box_radius = data.get(valid_flags["box_radius"], None)
        self.box = data.get(valid_flags["box"], None)
        self.native = data.get(valid_flags["native"], "")
        self.atom_dist = data.get(valid_flags["atom_dist"], None)
        self.debug = data.get(valid_flags["debug"], None)
        self.folder = data.get(valid_flags["folder"], None)
        self.output = data.get(valid_flags["output"], None)
        self.randomize = data.get(valid_flags["randomize"], None)
        self.full = data.get(valid_flags["full"], None)
        self.proximityDetection = data.get(valid_flags["proximityDetection"], None)
        self.poses = data.get(valid_flags["poses"], None)
        self.precision_glide = data.get(valid_flags["precision_glide"], None)
        self.msm = data.get(valid_flags["msm"], None)
        self.precision = data.get(valid_flags["precision"], None)
        self.clust = data.get(valid_flags["clust"], None)
        self.restart = data.get(valid_flags["restart"], None)
        self.lagtime = data.get(valid_flags["lagtime"], None)
        self.msm_clust = data.get(valid_flags["msm_clust"], None)
        self.rescoring = data.get(valid_flags["rescoring"], None)
        self.in_out = data.get(valid_flags["in_out"], None)
        self.in_out_soft = data.get(valid_flags["in_out_soft"], None)
        self.exit = data.get(valid_flags["exit"], None)
        self.exit_value = data.get(valid_flags["exit_value"], None)
        self.exit_condition = data.get(valid_flags["exit_condition"], None)
        self.exit_trajnum = data.get(valid_flags["exit_trajnum"], None)
        self.waters = data.get(valid_flags["waters"], None)
        self.water_freq = data.get(valid_flags["water_freq"], None)
        self.water_center = data.get(valid_flags["water_center"], None)
        self.water_temp = data.get(valid_flags["water_temp"], None)
        self.water_overlap = data.get(valid_flags["water_overlap"], None)
        self.water_constr = data.get(valid_flags["water_constr"], None)
        self.water_trials = data.get(valid_flags["water_trials"], None)
        self.water_radius = data.get(valid_flags["water_radius"], None)
        self.induced_fit_exhaustive = data.get(
            valid_flags["induced_fit_exhaustive"], None
        )
        self.induced_fit_fast = data.get(valid_flags["induced_fit_fast"], None)
        self.frag = data.get(valid_flags["frag"], None)
        self.ca_constr = data.get(valid_flags["ca_constr"], None)
        self.ca_interval = data.get(valid_flags["ca_interval"], None)
        self.constraint_level = data.get(valid_flags["constraint_level"], None)
        self.terminal_constr = data.get(valid_flags["terminal_constr"], None)
        self.one_exit = data.get(valid_flags["one_exit"], None)
        self.box_type = data.get(valid_flags["box_type"], None)
        self.box_metric = data.get(valid_flags["box_metric"], None)
        self.time = data.get(valid_flags["time"], None)
        self.nosasa = data.get(valid_flags["nosasa"], None)
        self.sasa = data.get(valid_flags["sasa"], None)
        self.perc_sasa = data.get(valid_flags["perc_sasa"], None)
        self.seed = data.get(valid_flags["seed"], None)
        self.pdb = data.get(valid_flags["pdb"], None)
        self.log = data.get(valid_flags["log"], None)
        self.nonrenum = data.get(valid_flags["nonrenum"], None)
        self.pele_exec = data.get(valid_flags["pele_exec"], None)
        self.pele_data = data.get(valid_flags["pele_data"], None)
        self.pele_documents = data.get(valid_flags["pele_documents"], None)
        self.pca = data.get(valid_flags["pca"], None)
        self.anm_direction = data.get(valid_flags["anm_direction"], None)
        self.anm_mix_modes = data.get(valid_flags["anm_mix_modes"], None)
        self.anm_picking_mode = data.get(valid_flags["anm_picking_mode"], None)
        self.anm_displacement = data.get(valid_flags["anm_displacement"], None)
        self.anm_modes_change = data.get(valid_flags["anm_modes_change"], None)
        self.anm_num_of_modes = data.get(valid_flags["anm_num_of_modes"], None)
        self.anm_relaxation_constr = data.get(
            valid_flags["anm_relaxation_constr"], None
        )
        self.remove_constraints = data.get(valid_flags["remove_constraints"], None)
        self.pca_traj = data.get(valid_flags["pca_traj"], None)
        self.perturbation = data.get(valid_flags["perturbation"], None)
        self.binding_energy = data.get(valid_flags["binding_energy"], None)
        self.parameters = data.get(valid_flags["parameters"], None)
        self.analyse = data.get(valid_flags["analyse"], None)
        self.selection_to_perturb = data.get(valid_flags["selection_to_perturb"], None)
        self.mae = data.get(valid_flags["mae"], None)
        self.constrain_core = data.get(valid_flags["constrain_core"], None)
        self.constrain_core_spring = data.get(
            valid_flags["constrain_core_spring"], 50.0
        )
        self.skip_ligand_prep = data.get(valid_flags["skip_ligand_prep"], None)
        self.spawning_condition = data.get(valid_flags["spawning_condition"], None)
        self.external_constraints = data.get(valid_flags["external_constraints"], [])
        self.only_analysis = data.get(valid_flags["only_analysis"], False)
        self.overwrite = data.get(valid_flags["overwrite"], True)
        self.analysis_nclust = data.get(valid_flags["analysis_nclust"], 10)
        self.te_column = data.get(valid_flags["te_column"], 4)
        self.be_column = data.get(valid_flags["be_column"], 5)
        self.limit_column = data.get(valid_flags["limit_column"], 6)
        self.com = data.get(valid_flags["com"], None)
        self.pele_license = data.get(valid_flags["pele_license"], None)
        self.schrodinger = data.get(valid_flags["schrodinger"], None)
        self.no_check = data.get(valid_flags["no_check"], False)
        self.cleanup = data.get(valid_flags["cleanup"], False)
        self.water_empty_selector = data.get(valid_flags["water_empty_selector"], False)
        self.polarize_metals = data.get(valid_flags["polarize_metals"], False)
        self.polarization_factor = data.get(valid_flags["polarization_factor"], 2)
        self.interaction_restrictions = data.get(valid_flags["interaction_restrictions"], None)
        self.inter_step_logger = data.get(valid_flags["inter_step_logger"], None)
        self.singularity_exec = data.get(valid_flags["singularity_exec"], None)

        # Metal constraints
        self.permissive_metal_constr = data.get(
            valid_flags["permissive_metal_constr"], False
        )
        self.constrain_all_metals = data.get(valid_flags["constrain_all_metals"], False)
        self.no_metal_constraints = data.get(valid_flags["no_metal_constraints"], False)

        # Frag
        self.frag_run = data.get(valid_flags["frag_run"], True)
        self.frag_core = data.get(valid_flags["frag_core"], False)
        self.frag_input = data.get(valid_flags["frag_input"], False)
        self.frag_ligands = data.get(valid_flags["frag_ligands"], False)
        self.growing_steps = data.get(valid_flags["growing_steps"], False)
        self.frag_steps = data.get(valid_flags["frag_steps"], False)
        self.frag_eq_steps = data.get(valid_flags["frag_eq_steps"], False)
        self.protocol = data.get(valid_flags["protocol"], None)
        self.frag_ai = data.get(valid_flags["frag_ai"], False)
        self.frag_ai_iterations = data.get(valid_flags["frag_ai_iterations"], False)
        self.chain_core = data.get(valid_flags["chain_core"], False)
        self.frag_restart = data.get(valid_flags["frag_restart"], False)
        self.frag_criteria = data.get(valid_flags["frag_criteria"], False)
        self.frag_output_folder = data.get(valid_flags["frag_output_folder"], False)
        self.frag_cluster_folder = data.get(valid_flags["frag_cluster_folder"], False)
        self.frag_library = data.get(valid_flags["frag_library"], None)
        self.frag_core_atom = data.get(valid_flags["frag_core_atom"], None)
        self.analysis_to_point = data.get(valid_flags["analysis_to_point"], None)
        self.fragment_atom = data.get(valid_flags["fragment_atom"], None)

        # PPI
        self.n_components = data.get(valid_flags["n_components"], None)
        self.ppi = data.get(valid_flags["ppi"], None)
        self.center_of_interface = data.get(valid_flags["center_of_interface"], None)
        self.protein = data.get(valid_flags["protein"], None)
        self.ligand_pdb = data.get(valid_flags["ligand_pdb"], None)
        self.skip_refinement = data.get(valid_flags["skip_refinement"], None)
        self.n_waters = data.get(valid_flags["n_waters"], 0)

        # site_finder
        self.site_finder = data.get(valid_flags["site_finder"], None)
        self.skip_refinement = data.get(valid_flags["skip_refinement"], None)
        self.site_finder_local = data.get(valid_flags["site_finder_local"], None)
        self.site_finder_global = data.get(valid_flags["site_finder_global"], None)

        # RNA
        self.rna = data.get(valid_flags["rna"], None)

        # GPCR
        self.gpcr_orth = data.get(valid_flags["gpcr_orth"], None)
        self.orthosteric_site = data.get(valid_flags["orthosteric_site"], None)
        self.initial_site = data.get(valid_flags["initial_site"], None)

        # OUTIN
        self.final_site = data.get(valid_flags["final_site"], None)

        # Mutagenesis
        self.saturated_mutagenesis = data.get(valid_flags["saturated_mutagenesis"], None)
        self.cpus_per_mutation = data.get(valid_flags["cpus_per_mutation"], None)

        # Analysis
        self.clustering_method = data.get(valid_flags["clustering_method"], None)
        self.bandwidth = data.get(valid_flags["bandwidth"], None)
        self.kde = data.get(valid_flags["kde"], None)
        self.kde_structs = data.get(valid_flags["kde_structs"], None)
        self.max_top_clusters = data.get(valid_flags["max_top_clusters"], None)
        self.min_population = data.get(valid_flags["min_population"], None)
        self.max_top_poses = data.get(valid_flags["max_top_poses"], None)
        self.top_clusters_criterion = data.get(valid_flags["top_clusters_criterion"], None)
        self.cluster_representatives_criterion = data.get(valid_flags["cluster_representatives_criterion"], None)

        if self.test:
            warnings.warn(
                "WARNING: This simulation is a test do not use the input files to run production simulations"
            )
            self.cpus = 5
            self.pele_steps = self.steps = 1
            self.iterations = 1
            self.min_freq = 0
            self.anm_freq = 0
            self.sidechain_freq = 0
            self.n_components = 3
            self.temperature = self.temp = 10000
            self.n_components = 3
            self.analysis_nclust = 4
            self.max_top_clusters = 4
            self.cpus_per_mutation = 2


@dataclass
class Most_Similar_Flag:
    name: str

    def calculate_distance(self, key):
        self.distance = SequenceMatcher(None, self.name, key).ratio()
