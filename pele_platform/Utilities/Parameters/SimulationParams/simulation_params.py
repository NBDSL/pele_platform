import random
import os
import glob

import pele_platform.constants.constants as cs
from pele_platform.Utilities.Parameters.SimulationParams.MSMParams import msm_params
from pele_platform.Utilities.Parameters.SimulationParams.GlideParams import glide_params
from pele_platform.Utilities.Parameters.SimulationParams.BiasParams import bias_params
from pele_platform.Utilities.Parameters.SimulationParams.InOutParams import inout_params
from pele_platform.Adaptive.interaction_restrictions import InteractionRestrictionsBuilder
from pele_platform.Utilities.Parameters.SimulationParams.PCA import pca
from pele_platform.Utilities.Parameters.SimulationParams.site_finder import site_finder
from pele_platform.Utilities.Parameters.SimulationParams.PPI import ppi
import pele_platform.Utilities.Helpers.helpers as hp

LOGFILE = '"simulationLogPath" : "$OUTPUT_PATH/logFile.txt",'


class SimulationParams(
    msm_params.MSMParams,
    glide_params.GlideParams,
    bias_params.BiasParams,
    inout_params.InOutParams,
    pca.PCAParams,
    site_finder.SiteFinderParams,
    ppi.PPIParams
):

    def __init__(self, args):
        self.simulation_type(args)
        self.main_pele_params(args)
        self.singularity_params(args)
        self.main_adaptive_params(args)
        self.optative_params(args)
        self.system_preparation_params(args)
        self.ligand_params(args)
        self.anm_params(args)
        self.water_params(args)
        self.box_params(args)
        self.output_params(args)
        self.analysis_params(args)
        self.constraints_params(args)
        self.interaction_restrictions_params(args)

        # Create all simulation types (could be more efficient --> chnage in future)
        super().generate_msm_params(args)
        super().generate_glide_params(args)
        super().generate_bias_params(args)
        super().generate_inout_params(args)
        super().generate_pca_params(args)
        super().generate_site_finder_params(args)
        super().generate_ppi_params(args)
        # rna.RNAParams.__init__(self, args)

    def simulation_type(self, args):
        self.adaptive = True if args.package in ["site_finder", "adaptive", "PPI"] else None
        self.frag_pele = True if args.package == "frag" else None
        # Trick to let frag handle control fodler parameters --> Improve
        self.complexes = "$PDB" if self.frag_pele else "$COMPLEXES"
        self.frag_pele_steps = "$STEPS" if self.frag_pele else "$PELE_STEPS"
        self.output_path = "$RESULTS_PATH" if self.frag_pele else "$OUTPUT_PATH"

    def main_pele_params(self, args):
        if "*" in args.system:
            self.system = glob.glob(args.system)[0]
            self.input = glob.glob(args.system) if not args.input else args.input
        else:
            self.system = args.system
            self.input = args.input
        self.residue = args.residue
        self.chain = args.chain
        if self.adaptive:
            assert (
                    self.system and self.residue and self.chain
            ), "User must define input, residue and chain"
        self.debug = args.debug if args.debug else False
        self.pele_steps = (
            args.pele_steps
            if args.pele_steps
            else self.simulation_params.get("pele_steps", 8)
        )
        self.logfile = LOGFILE if args.log else ""
        self.license = (
            args.pele_license
            if args.pele_license
            else cs.DEFAULT_PELE_LICENSE
        )
        self.anm_freq = (
            args.anm_freq
            if args.anm_freq is not None
            else self.simulation_params.get("anm_freq", 4)
        )
        self.sidechain_freq = (
            args.sidechain_freq
            if args.sidechain_freq is not None
            else self.simulation_params.get("sidechain_freq", 2)
        )
        self.min_freq = (
            args.min_freq
            if args.min_freq is not None
            else self.simulation_params.get("min_freq", 1)
        )
        self.water_freq = (
            args.water_freq
            if args.water_freq is not None
            else self.simulation_params.get("water_freq", 1)
        )
        self.temperature = (
            args.temperature
            if args.temperature
            else self.simulation_params.get("temperature", 1500)
        )
        self.sidechain_resolution = (
            args.sidechain_resolution
            if args.sidechain_resolution
            else self.simulation_params.get("sidechain_resolution", 30)
        )
        self.proximityDetection = (
            "false"
            if args.proximityDetection is False
            else self.simulation_params.get("proximityDetection", "true")
        )
        self.steric_trials = (
            args.steric_trials
            if args.steric_trials
            else self.simulation_params.get("steric_trials", 250)
        )
        self.overlap_factor = (
            args.overlap_factor
            if args.overlap_factor
            else self.simulation_params.get("overlap_factor", 0.65)
        )
        self.steering = (
            args.steering
            if args.steering
            else self.simulation_params.get("steering", 0)
        )
        self.perturbation = (
            ""
            if args.perturbation is False
            else self.simulation_params.get("perturbation", cs.PERTURBATION)
        )
        self.perturbation_params(args)
        self.com = (
            args.com
            if args.com
            else self.simulation_params.get("COMligandConstraint", 0)
        )

    def anm_params(self, args):
        self.anm_displacement = (
            args.anm_displacement
            if args.anm_displacement
            else self.simulation_params.get("anm_displacement", 0.75)
        )
        self.anm_modes_change = (
            args.anm_modes_change
            if args.anm_modes_change
            else self.simulation_params.get("anm_modes_change", 4)
        )
        self.anm_direction = (
            args.anm_direction
            if args.anm_direction
            else self.simulation_params.get("anm_direction", "random")
        )
        self.anm_mix_modes = (
            args.anm_mix_modes
            if args.anm_mix_modes
            else self.simulation_params.get(
                "anm_mix_modes", "mixMainModeWithOthersModes"
            )
        )
        self.anm_picking_mode = (
            args.anm_picking_mode
            if args.anm_picking_mode
            else self.simulation_params.get("anm_picking_mode", "RANDOM_MODE")
        )
        self.anm_num_of_modes = (
            args.anm_num_of_modes
            if args.anm_num_of_modes
            else self.simulation_params.get("anm_num_of_modes", 6)
        )
        self.anm_relaxation_constr = (
            args.anm_relaxation_constr
            if args.anm_relaxation_constr
            else self.simulation_params.get("anm_relaxation_constr", 0.5)
        )
        self.remove_constraints = (
            args.remove_constraints
            if args.remove_constraints is not None
            else self.simulation_params.get("remove_constraints", False)
        )

    def perturbation_params(self, args):
        if self.perturbation:
            self.selection_to_perturb = (
                args.selection_to_perturb
                if args.selection_to_perturb
                else self.simulation_params.get(
                    "selection_to_perturb", cs.SELECTION_TO_PERTURB
                )
            )
            self.parameters = (
                args.parameters
                if args.parameters
                else self.simulation_params.get("params", True)
            )
            self.ligand = cs.LIGAND if self.perturbation else ""
            self.binding_energy = (
                args.binding_energy
                if args.binding_energy
                else self.simulation_params.get("binding_energy", cs.BE)
            )
            self.sasa = (
                args.sasa if args.sasa else self.simulation_params.get("sasa", cs.SASA)
            )
        else:
            self.selection_to_perturb = ""
            self.parameters = ""
            self.ligand = ""
            self.binding_energy = ""
            self.sasa = ""

    def main_adaptive_params(self, args):
        self.spawning = (
            args.spawning
            if args.spawning
            else self.simulation_params.get("spawning_type", "independent")
        )
        self.spawning_condition = (
            args.spawning_condition
            if args.spawning_condition
            else self.simulation_params.get("spawning_condition", None)
        )
        self.spawning_condition = (
            '"condition": "{}",'.format(self.spawning_condition)
            if self.spawning_condition
            else ""
        )
        self.density = (
            args.density
            if args.density
            else self.simulation_params.get("density", "null")
        )
        self.simulation_type = (
            args.simulation_type
            if args.simulation_type
            else self.simulation_params.get("simulation_type", "pele")
        )
        iterations = 1 if self.spawning == "independent" else 30
        self.iterations = (
            args.iterations
            if args.iterations
            else self.simulation_params.get("iterations", iterations)
        )
        self.clust_type = (
            args.clust_type
            if args.clust_type
            else self.simulation_params.get("clust_type", "rmsd")
        )
        self.cluster_values = (
            args.cluster_values
            if args.cluster_values
            else self.simulation_params.get("cluster_values", "[1.75, 2.5, 4, 6]")
        )
        self.cluster_conditions = (
            args.cluster_conditions
            if args.cluster_conditions
            else self.simulation_params.get("cluster_conditions", "[1, 0.6, 0.4, 0.0]")
        )
        self.seed = args.seed if args.seed else random.randrange(1, 70000)
        self.templates = os.path.abspath(
            os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates")
        )
        self.usesrun = "true" if args.usesrun else "false"
        mpi_params_name = "srunParameters" if args.usesrun else "mpiParameters"
        self.mpi_params = (
            f'"{mpi_params_name}": "{args.mpi_params}",' if args.mpi_params else ""
        )

    def optative_params(self, args):
        self.forcefield = (
            args.forcefield
            if args.forcefield
            else self.simulation_params.get("forcefield", "OPLS2005")
        )
        self.solvent = (
            args.solvent
            if args.solvent
            else self.simulation_params.get("solvent", "VDGBNP")
        )
        self.verbose = (
            "true" if args.verbose else self.simulation_params.get("verbose", "false")
        )
        self.cpus = args.cpus = (
            args.cpus if args.cpus else self.simulation_params.get("cpus", 60)
        )

        self.cpus_per_mutation = args.cpus_per_mutation if args.cpus_per_mutation else 0

        self.restart = (
            args.restart
            if args.restart
            else self.simulation_params.get("restart", "all")
        )
        self.test = args.test
        # +1 to avoid being 0
        self.equil_steps = (
            int(args.eq_steps / self.cpus) + 1
            if args.eq_steps
            else self.simulation_params.get("equilibration_steps", 1)
        )
        self.equilibration = (
            "true"
            if args.equilibration
            else self.simulation_params.get("equilibration", "false")
        )
        self.equilibration_mode = (
                args.equilibration_mode 
                if args.equilibration_mode
                else self.simulation_params.get("equilibration_mode", "equilibrationSelect")
        )
        self.adaptive_restart = args.adaptive_restart
        self.poses = (
            args.poses
            if args.poses
            else self.simulation_params.get("poses", self.cpus - 1)
        )
        self.pele_exec = (
            args.pele_exec if args.pele_exec else os.path.join(cs.PELE, "bin/Pele_mpi")
        )
        self.pele_data = (
            args.pele_data if args.pele_data else os.path.join(cs.PELE, "Data")
        )
        self.pele_documents = (
            args.pele_documents
            if args.pele_documents
            else os.path.join(cs.PELE, "Documents")
        )
        self.polarize_metals = args.polarize_metals if args.polarize_metals else False
        self.polarization_factor = (
            args.polarization_factor if args.polarization_factor else 2.0
        )
        self.skip_refinement = args.skip_refinement if args.skip_refinement else False
        self.bandwidth = args.bandwidth if args.bandwidth else self.simulation_params.get(
            "bandwidth", 2.5)
        self.clustering_method = args.clustering_method if args.clustering_method else self.simulation_params.get(
            "clustering_method", "meanshift")

    def system_preparation_params(self, args):
        self.skip_prep = (
            args.skip_prep
            if args.skip_prep
            else self.simulation_params.get("skip_prep", False)
        )
        self.nonstandard = (
            args.nonstandard
            if args.nonstandard
            else self.simulation_params.get("nonstandard", [])
        )
        self.constraints = None
        self.external_constraints = (
            hp.retrieve_constraints_for_pele(args.external_constraints, self.system)
            if args.external_constraints
            else []
        )
        self.permissive_metal_constr = (
            args.permissive_metal_constr if args.permissive_metal_constr else []
        )
        self.constrain_core = (
            args.constrain_core
            if args.constrain_core
            else self.simulation_params.get("constrain_core", None)
        )
        self.constrain_core_spring = (
            args.constrain_core_spring if args.constrain_core_spring else 50.0
        )
        self.no_ppp = (
            args.no_ppp if args.no_ppp else self.simulation_params.get("no_ppp", False)
        )

    def ligand_params(self, args):
        self.spython = os.path.join(cs.SCHRODINGER, "utilities/python")
        if not os.path.exists(self.spython):
            self.spython = os.path.join(cs.SCHRODINGER, "run")
        self.mae_lig = (
            args.mae_lig
            if args.mae_lig
            else self.simulation_params.get("mae_lig", None)
        )
        self.external_template = (
            args.template
            if args.template
            else self.simulation_params.get("template", [])
        )
        self.external_rotamers = (
            args.rotamers
            if args.rotamers
            else self.simulation_params.get("rotamers", [])
        )
        self.skip_ligand_prep = (
            args.skip_ligand_prep
            if args.skip_ligand_prep
            else self.simulation_params.get("args.skip_ligand_prep", [])
        )
        self.core = args.core
        self.n = args.n
        self.mtor = args.mtor
        self.forcefield = args.forcefield
        self.mae_lig = args.mae_lig
        self.lig = self.mae_lig if self.mae_lig else "{}.mae".format(self.residue)
        self.gridres = args.gridres

    def water_params(self, args):
        self.water_temp = (
            args.water_temp
            if args.water_temp
            else self.simulation_params.get("water_temp", 5000)
        )
        self.water_overlap = (
            args.water_overlap
            if args.water_overlap
            else self.simulation_params.get("water_overlap", 0.78)
        )
        self.water_constr = (
            args.water_constr
            if args.water_constr
            else self.simulation_params.get("water_constr", 0)
        )
        self.water_trials = (
            args.water_trials
            if args.water_trials
            else self.simulation_params.get("water_trials", 10000)
        )

        self.allow_empty_selectors = (
            '"allowEmptyWaterSelectors": true,' if args.water_empty_selector else ""
        )
        self.water_arg = (
            hp.retrieve_all_waters(self.system)
            if args.waters == "all_waters"
            else args.waters
        )  # IDs of waters
        self.n_waters = (
            args.n_waters
            if args.n_waters
            else self.simulation_params.get("n_waters", 0)
        )
        self.waters = (
            args.waters if args.waters else self.simulation_params.get("waters", [])
        )
        self.water_radius = args.water_radius if args.water_radius else 6
        self.water_center = None
        self.water = ""
        self.water_energy = None

    def box_params(self, args):
        self.box_radius = (
            args.box_radius
            if args.box_radius
            else self.simulation_params.get("box_radius", None)
        )
        if args.box_center:
            if ":" in args.box_center:
                self.box_center = [
                    str(coord)
                    for coord in hp.get_coords_from_residue(
                        self.system, args.box_center
                    )
                ]
                self.box_center = "[" + ", ".join(self.box_center) + "]"
            else:
                self.box_center = [str(x) for x in args.box_center]
                self.box_center = (
                        "[" + ",".join([str(coord) for coord in self.box_center]) + "]"
                )
        else:
            self.box_center = self.simulation_params.get("box_center", None)

    def output_params(self, args):
        self.folder = args.folder
        self.output = args.output if args.output is not None else "output"
        self.report_name = args.report_name if args.report_name else "report"
        self.traj_name = args.traj_name if args.traj_name else "trajectory.pdb"
        self.xtc = self.traj_name.endswith(".xtc")
        self.pdb = self.traj_name.endswith(".pdb")
        self.inter_step_logger = cs.INTERSTEPLOGGER if args.inter_step_logger else ""

    def analysis_params(self, args):
        self.analyse = args.analyse if args.analyse is not None else True
        self.mae = args.mae if args.mae else False
        self.only_analysis = args.only_analysis
        self.analysis_nclust = args.analysis_nclust
        self.overwrite = args.overwrite
        self.be_column = args.be_column
        self.te_column = args.te_column
        self.limit_column = args.limit_column
        self.kde = args.kde if args.kde is not None else False
        self.kde_structs = args.kde_structs if args.kde_structs else 1000
        self.min_population = args.min_population if args.min_population is not None else 0.01
        self.max_top_clusters = args.max_top_clusters if args.max_top_clusters is not None else 8
        self.max_top_poses = args.max_top_poses if args.max_top_poses is not None else 100
        self.top_clusters_criterion = args.top_clusters_criterion if args.top_clusters_criterion is not None else "interaction_25_percentile"
        self.cluster_representatives_criterion = args.cluster_representatives_criterion if args.cluster_representatives_criterion is not None else "interaction_5_percentile"

    def constraints_params(self, args):
        """
        Sets backbone constraints parameters, such as interval and spring constants. It prioritises the user-defined
        values, then checks default values for specified constraints level (if any), then checks simulation_params for
        a specific package to finally fall back to level 1 default.
        Args:
            args (EnviroBuilder): arguments defined by the user in input.yaml

        Returns:
            Backbone CA and terminal CA spring constants, backbone CA interval.
        """
        constraints_flags = [
            "ca_interval",
            "ca_constr",
            "terminal_constr",
        ]  # list of relevant flags to check

        # Checks, if the user specified one of the built-in constraint levels.
        self.constraint_level = (
            args.constraint_level if args.constraint_level is not None else None
        )
        self.level_params = cs.constraint_levels.get(self.constraint_level, {})

        # Dictionary with default constraints params, i.e. level 1
        defaults = cs.constraint_levels.get(1, {})

        for flag in constraints_flags:
            flag_value = getattr(args, flag, None)
            if not flag_value:
                flag_value = self.level_params.get(
                    flag, self.simulation_params.get(flag, defaults[flag])
                )
            setattr(self, flag, flag_value)

    def interaction_restrictions_params(self, args):
        """
        Sets parameters for interaction restrictions.
        Fills the pele_params INTERACTION_RESTRICTIONS template with an additional section of parameters change.
        """
        if args.interaction_restrictions:
            restrictions = InteractionRestrictionsBuilder()
            restrictions.parse_interaction_restrictions(self.system, args.interaction_restrictions)
            self.met_interaction_restrictions = restrictions.metrics_to_json()
            self.interaction_restrictions = restrictions.conditions_to_json()
            self.parameters = restrictions.fill_template(self.parameters)
        else:
            self.met_interaction_restrictions = ""
            self.interaction_restrictions = ""

    def singularity_params(self, args):
        """
        Sets parameters for singularity containers.
        """
        args.mpi_params = args.singularity_exec if args.singularity_exec else args.mpi_params
        if args.singularity_exec:
            args.pele_exec = "Pele_mpi" if not self.frag_pele else args.mpi_params + " Pele_mpi"
