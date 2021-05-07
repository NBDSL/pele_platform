from dataclasses import dataclass
import os
import AdaptivePELE.adaptiveSampling as ad
from pele_platform.Utilities.Helpers import helpers, template_builder
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Utilities.Helpers.water as wt


@dataclass
class SimulationBuilder(template_builder.TemplateBuilder):

    adaptive_file: str
    pele_file: str
    topology: str

    def generate_inputs(self, parameters, water_obj):
        """
        It generates the input files for Adaptive PELE, according to the
        parameters previously generated.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        water_obj : a WaterIncluder object
            The parameters for aquaPELE, if applicable
        """
        # Fill in simulation template with 
        # simulation parameters specified in  

        # Fill two times because we have flags inside flags
        self.fill_pele_template(parameters, water_obj)
        self.fill_pele_template(parameters, water_obj)
        self.fill_adaptive_template(parameters)
        self.fill_adaptive_template(parameters)

    def fill_pele_template(self, env: pv.ParametersBuilder, water_obj: wt.WaterIncluder) -> None:
        # Fill in PELE template
        self.pele_keywords = { "PERTURBATION": env.perturbation, "SELECTION_TO_PERTURB": env.selection_to_perturb,
                        "BE": env.binding_energy, "SASA": env.sasa,
                        "LOGFILE": env.logfile, "NATIVE": env.native, "FORCEFIELD": env.forcefield, "CHAIN": env.chain, 
                        "CONSTRAINTS": "\n".join(env.constraints), "CPUS":env.cpus,
                        "LICENSES": env.license, "BOX_RADIUS": env.box_radius, "BOX_CENTER": env.box_center,
                        "SASA_min": env.sasa_min, "SASA_max": env.sasa_max,
                        "WATER_RADIUS": water_obj.water_radius, "WATER_CENTER": water_obj.water_center, "WATER": water_obj.water_line,
                        "WATER_ENERGY": water_obj.water_energy, "METRICS": env.metrics, "REPORT_NAME": env.report_name, 
                        "TRAJ_NAME": env.traj_name,
                        "SOLVENT": env.solvent, "PARAMETERS": env.parameters, "SIDECHAIN_RESOLUTION": env.sidechain_resolution,
                        "OVERLAP": env.overlap_factor, "STERIC_TRIALS": env.steric_trials, "TEMPERATURE": env.temperature, 
                        "MIN_FREQ": env.min_freq, "SIDECHAIN_FREQ": env.sidechain_freq, "ANM_FREQ": env.anm_freq, "BOX" : env.box, "PROXIMITY": env.proximityDetection, "WATER_FREQ": env.water_freq, "VERBOSE": env.verbose, "ANM_DISPLACEMENT": env.anm_displacement, "ANM_MODES_CHANGE": env.anm_modes_change, "ANM_DIRECTION": env.anm_direction, "ANM_MIX_MODES": env.anm_mix_modes, "ANM_PICKING_MODE": env.anm_picking_mode,
                        "ANM_NUM_OF_MODES": env.anm_num_of_modes, "ANM_RELAXATION_CONST": env.anm_relaxation_constr,
                        "PCA": env.pca, "COMPLEXES": env.complexes, "PELE_STEPS": env.frag_pele_steps, "OUTPUT_PATH": env.output_path,
                        "COM": env.com, "STEERING": env.steering, "MET_INTERACTION_RESTRICTIONS":env.met_interaction_restrictions,
                        "INTERACTION_RESTRICTIONS": env.interaction_restrictions, "INTER_STEP_LOGGER": env.inter_step_logger}

        super(SimulationBuilder, self).__init__(self.pele_file, self.pele_keywords)

    def fill_adaptive_template(self, env: pv.ParametersBuilder) -> None:
        # Fill in adaptive template
        self.adaptive_keywords = { "RESTART": env.adaptive_restart, "OUTPUT": env.output, "INPUT":env.adap_ex_input,
                "CPUS":env.cpus, "PELE_CFILE": os.path.basename(self.pele_file), "LIG_RES": env.residue, "SEED": env.seed, "EQ_STEPS": env.equil_steps,
                "EQUILIBRATION":env.equilibration, "EQUILIBRATION_MODE": env.equilibration_mode,"EPSILON": env.epsilon, 
                "BIAS_COLUMN": env.bias_column, "ITERATIONS": env.iterations, 
                "PELE_STEPS": env.pele_steps, "REPORT_NAME": env.report_name, "SPAWNING_TYPE": env.spawning, "DENSITY": env.density,
                "SIMULATION_TYPE": env.simulation_type, "CLUSTER_VALUES": env.cluster_values, "CLUSTER_CONDITION": env.cluster_conditions,
                "UNBINDING": env.unbinding_block, "USESRUN": env.usesrun, "LIGAND": env.ligand,
                "PELE_BIN": env.pele_exec, "PELE_DATA": env.pele_data, "PELE_DOCUMENTS": env.pele_documents,
                "CONDITION": env.spawning_condition, "MPI_PARAMS": env.mpi_params, "CLUST_TYPE": env.clust_type}
        super(SimulationBuilder, self).__init__(self.adaptive_file, self.adaptive_keywords)


    def run(self, hook=False) -> None:
        # Launch montecarlo simulation
        with helpers.cd(os.path.dirname(self.adaptive_file)):
                ad.main(self.adaptive_file)

