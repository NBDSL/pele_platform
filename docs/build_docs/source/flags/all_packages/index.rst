All packages
===============

General settings
---------------------

Configure the settings of the simulation and the path to all dependencies in case of need (non-default installation).

- **test**: Run a quick test to check the simulation works (~2 min). **Never use the control files from the test as input for a production simulation as temperature, ANM and minimization are twicked to make simulation faster!!!!**

- **usesrun**: Use srun binary to run PELE. Only when using intel processors.

- **pele_exec**: Use a pele executable that is not the default one. **Needs to be used with pele_data and pele_documents**. default: $PELE/bin/Pele_mpi

- **pele_data**: Use a pele data folder that is not the default one. default: $PELE/Data

- **pele_documents**: Use a pele documents folder that is not the default one. default: $PELE/Documents

- **singularity_exec**: Use a singularity container that contains the pele executable. **Cannot be used with pele_exec**.

- **pele_license**: Use a pele_license path that is not the default one. default: $PELE/licenses

- **schrodinger**: Use a schrodinger path that is not the default one. default: $SCHRODINGER

..  code-block:: yaml


  test: true
  usesrun: false
  pele_exec: "/home/pele/bin/Pele_mpi"
  pele_data: "/home/pele/Data/"
  pele_documents: "/home/pele/Documents/"
  singularity_exec: "/home/pele/pele_release.sif"
  pele_license: "/home/pele/licenses"
  schrodinger: "/home/pele/schrodinger2020-1/"


Simulation parameters
----------------------

- **seed**: Seed of the job for reproducibility. Default=12345

- **log**: Retrieve PELE logfiles during simulation. Default=False

- **verbose**: Set to true to activate verbose mode in PELE. DEfault=False

- **anm_freq**: Every how many steps to perform anm. Default=4

- **anm_displacement**: Angstrom to displace carbon alphas in each ANM movement. Default=0.75

- **anm_modes_change**: Number of steps before we change to a new normal mode movement. Default=4

- **sidechain_freq**: Every how many steps to perform sidechain sampling. Default=2

- **min_freq**: Every how many steps to perform minimization. Default=1

- **water_freq**: Every how many steps to perform water perturbation. Default=1

- **temperature**: Temperature of the simulation. Default=1500

- **solvent**: Solvent of the simulation. (OBC or VDGBNP). Please, take into account that any OpenFF force field can only run with the OBC solvent model. Default=VDGBNP (unless for the OpenFF, which the default is OBC)

- **sidechain_res**: Receptor sidechain resolution. Default=10

- **overlap_factor**: Vanderwals overlap factor (More in PELE docs). Default=0.65

- **steric_trials**: Number of steric trials (More in PELE docs). Default=250

..  code-block:: yaml

  seed: 312312
  log: true
  verbose: true
  anm_freq: 4
  anm_displacement: 0.5
  anm_modes_change: 3
  sidechain_freq: 2
  min_freq: 1
  water_freq: 1
  temperature: 1500
  solvent: "VDGBNP"
  sidechain_res: 30
  overlap_factor: 0.65
  steric_trials: 250


PELE++ parameters
-------------------

**These flags are exclusive of the PELE modes not fragPELE**

- **iterations**: Adaptive epochs to run. Set to 1 by default if using PELE

- **steps**: Pele steps in each iteration

- **debug**: Use this flag to only create the inputs of the simulation. No simulation is run. (Usefull to transport it to another machine)

- **spawning**: Spawning type ([independent, inverselyProportional or epsilon so far]). Default: inverselyProportional

- **density**: Density type ([null, exitContinuous...]. More in AdaptivePELE docs). Default: null

- **cluster_values**: Clusterization values. More in AdaptivePELE. Default: Depending on simulation type

- **cluster_conditions**: Clusterization condition. More in AdaptivePELE. Default: Depending on simulation type

- **equilibration**: Whether to run initial equilibration or not. Default: false

- **equilibration_steps**: Equilibration steps. Default: 2

- **adaptive_restart**: Use adaptive restart with the working folder option to restart the simulation. Default: false

- **report**: Change the name of the report file. Default: report

- **traj**: Change the name of the trajectory file. Default: trajectory.pdb

..  code-block:: yaml

    iterations: 30
    steps: 12
    debug: true
    spawning: "epsilon"
    density: "exitContinuous"
    cluster_values: [2,3,4]
    cluster_conditions: [0.8, 0.6, 0.2]
    equilibration: false
    equilibration_steps: 10
    adaptive_restart: true
    working_folder: "folder_to_restart"
    report: report
    traj: trajectory.xtc

Receptor preparation
-----------------------

Configure the parameters of the PPP (Protein Pele Preparation)

- **skip_preprocess**: Skip protein pele preparation. Default: False

- **noTERs**: Don't include TERs on preparation. Used if PPP gets confuse with insertion codes or other. Default: False

- **charge_ters**: Charge terminals of the protein. Default: False

- **nonstandard**: List of names of nonstandard residues that will be omitted in protein pele preparation. Default=[]

- **prepwizard**: Run Prepwizard (Still on testing version). Default: False

..  code-block:: yaml

  preprocess_receptor: true
  noTERs: false
  charge_ters: false
  nonstandard:
    - TPO
  prepwizard: false


Ligand preparation
------------------

In order to run a simulation, PELE requires the following files for every non standard molecule (i.e. any non standard small molecule or residue):

    - **IMPACT template**: it contains the force field parameters. Please, check `this site <https://nostrumbiodiscovery.github.io/pele_docs/fileFormats.html#impact-template-file-format>`_ to get further information about it.
    - **rotamer library**: optional file containing the list of rotatable bonds to sample by the side chain perturbation algorithm. If missing, the flexibility of the corresponding molecule will not considered. Please, check `this site <https://nostrumbiodiscovery.github.io/pele_docs/fileFormats.html#sec-fileformats-ligandrotamers>`_
    - **solvent template**: some special solvents like "OBC" require extra parameters.

All these files can be automatically generated with `peleffy (PELE Force Field Yielder) <https://github.com/martimunicoy/peleffy>`_, one of the dependencies of the PELE Platform.
The following parameters below control the way how the PELE Platform will parametrize non standard molecules for you:

- **forcefield**: Forcefield used to parametrize hetero molecules, you can use one of:

        - "OPLS2005" (default)
        - "openff-1.3.0"
        - "openff-1.2.1"
        - "openff-1.2.0"
        - "openff-1.1.1"
        - "openff-1.1.0"
        - "openff-1.0.1"
        - "openff-1.0.0"

- **charge_parametrization_method**: The method to use to assign partial charges to atoms:

        - "gasteiger"
        - "am1bcc" (default when using any "OpenFF" force field)
        - "OPLS2005" (default when using "OPLS2005")

- **gridres**: Resolution of the rotamers when sampling them by the Side Chain prediction algorithm. Default=10 degrees

- **core**: List of PDB atom names that will be included as part of the rigid core. In case it is not specified, the algorithm will pick up a set of non-rotatable atoms centered in the molecular structure. Default=None

- **exclude_terminal_rotamers**: Exclude terminal rotamers during parametrization of non standard molecules if they belong to a small terminal group. Default=True

- **mae_lig**: External MAE file with quantum charges generated with Schrödinger suite. When supplied, any charge calculated internally in the platform will be replaced by the charges from this file. Default=None

..  code-block:: yaml

    charge_parametrization_method: "gasteiger"
    forcefield: "openff-1.3.0"
    gridres: 20
    core:
        - "O1"
        - "C1"
        - "C2"
        - "N1"

Alternatively, you can provide your own template and/or rotamer files as long as they follow PELE's naming convention
(see examples in the block code below). When templates

    - **templates**: External forcefield template files.

    - **rotamers**: External rotamer library files.

..  code-block:: yaml

  templates:
    - "/home/simulation_files/mgz"
    - "/home/simulation_files/ligz"
  rotamers:
    - "/home/simulation_files/MG.rot.assign"
    - "/home/simulation_files/LIG.rot.assign"

For more technical details about ligand parametrization, you can refer to the `PELE Force Field Yielder documentation <https://martimunicoy.github.io/peleffy/>`_.

Constraints
-----------

- **water_constr**: Water constraints. Default=5

- **constrain_smiles**: SMILES string to indicate what part of the molecule to constrain. Default=None

- **external_constraints**: You can specify 2 types of constraints: positional constraints or atom-atom constraints, e.g.

  - The positional constraints are given either by:
        - springConstant-atomnumber. i.e. "10-17"
        - springConstant-chain:resnum:atomname. i.e. "5-A:1:H"

  - The atom-atom constraints are specified either by:
        - springConstant-equilibriumDistance-atomnumber1-atomnumber2. i.e. "50-2.34-17-4159"
        - springConstant-equilibriumDistance-chain1:resnum1:atomname1-chain2:resnum2:atomname2. i.e. "50-2.34-A:1:H-L:1:C21"

- **remove_constraints**: Do not place constraints on the carbon-alpha of the protein. Default: False


..  code-block:: yaml

    water_constr: 5
    constrain_smiles: "C2CCC1CCCCC1C2"
    smiles_constr: 5
    external_constraints:
    - "10-17" #constrain of 10kcal/mol at atomnumber 17
    - "5-A:1:H" ##constrain of 10kcal/mol at atom with chain A residuenumber 1 and atomname H
    - "50-2.34-17-4159" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atomnumbers 17 & 4159
    - "50-2.34-A:1:H-L:1:C21" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    remove_constraints: true

Carbon-alpha constraints
+++++++++++++++++++++++++

Each package in the platform has its own predefined constraint parameters which are likely to be the best choice in each
type of study. However, the platform provides the users with several different levels of constraining the alpha carbons
of the protein backbone with varying spring constants and intervals:

- **level 0** - no constraints

- **level 1** - terminal CAs constrained with a spring constant of 5 kcal/mol, the rest of the CAs in the backbone with 0.5 kcal/mol at an interval of 10, i.e. every 10 residues (default)

- **level 2** - terminal CAs constrained at 5 kcal/mol, the rest of the CAs with 2.5 kcal/mol at the interval of 8 (default for the ``rescoring`` package)

- **level 3** - the whole backbone is constrained every 5 atoms with 5 kcal/mol (default for the ``gpcr_orth`` package)

We strongly suggest relying on the default settings for each package. However, in case of studying a system where the
defaults are not optimal (more flexibility or rigidity required), the users can change the level, for example:

..  code-block:: yaml

    constraint_level: 3

Alternatively, advanced users can manipulate the constraint parameters individually at their own risk, using the following flags:

- **terminal_constr** - sets the spring constant for the terminal C-alpha constraints, default = 5 kcal/mol

- **ca_constr** - sets the spring constant for the remaining C-alphas in the backbone, default = 0.5 kcal/mol

- **ca_interval** - interval at which the backbone C-alphas should be constrained, default = 10 (i.e. every 10 residues).

Take into account that specific modifiers of constraint parameters will prevail over the settings coming from the
constraints levels and those predefined in each package.

..  code-block:: yaml

    terminal_constr: 10.5
    ca_constr: 6.0
    ca_interval: 3

Metal constraints
+++++++++++++++++++++

Algorithm to automatically set metal constraints around the ligand.

- **no_metal_constraints**: Ignore all metals in the PDB file, no constraints will be set automatically. Default=False

- **permissive_metal_constr**: Expand the search for coordinated atoms by allowing 35% deviation from “ideal” angles. If the algorithm finds a valid geometry it will include the metal constraint into the simulation. Default=False

- **constrain_all_metals**: Constrain all atoms around the metal, regardless of the angles or coordination number. Default=False

- **external_constraints**: Set a manual constraint containing a metal atom to disable search for this particular metal. Default=[]


..  code-block:: yaml

    no_metal_constraints: true
    permissive_metal_constr: true
    constrain_all_metals: true
    external_constraints:
        - "50-2.34-A:1:H-L:1:MG" #constrain of 50kcal/mol with equilibrium distance of 2.34 between atoms with respective chain resnum and atomname
    constrain_core: "CN(C)C(=O)c1ccc(F)cc1"  # SMILES or SMARTS pattern
    constrain_core_spring: 30  # optional, default 50.0


Box parameters
---------------

Parameters to set the exploration Box:

- **box_radius**: Radius of the box. Default=[induced_fit (10), local_exploration (30), global_exploration (50)]

- **box_center**: Center of the box. Default=[indeuced_fit&local_exploration (CM of the ligand), global (calculater center)]


..  code-block:: yaml

  box_radius: 30
  box_center:
    - 20
    - 30
    - 50

Metal polarisation
-------------------

An optional flag to adjust charges on the metals by dividing them by certain factor.

- **polarize_metals** - adjust charges on the metals by dividing them by 2 (unless other value is set in polarization_factor)

- **polarization_factor** - factor by which the metal charges should be divided

..  code-block:: yaml

    polarize_metals: true
    polarization_factor: 2 # Mg2+ will have a charge of +1


Water perturbation
-----------------------

- **n_waters**: Number of waters to randomly add into your simulation and perturb. Default=0

- **waters**: Water molecules to be perturbed in AquaPELE steps. Users can indicate specific water IDs, e.g. "W:15" or select "all_waters" option to perturb all water molecules present in the system.

- **box_water**: Center of the box for the waters. Default: Centroid of the center of masses of all water molecules.

- **water_radius**: Radius of the water box. Default=7

- **water_trials**: Numerical trials on water perturbation. Default=10000

- **water_constr**: COM constrain applied to th water molecule after perturbation. Default=0

- **water_temp**: Temperature of the water perturbation step. Default=5000

- **water_overlap**: Overlap factor of water. Default=0.78


..  code-block:: yaml

    n_waters: 3 # Compulsory, if no water molecules are present in the system
    box_water:
    - 20
    - 30
    - 20
    water_radius: 8
    water_trials: 500
    water_constr: 0.5
    water_temp: 2000
    water_overlap: 0.5
    # waters: "all_waters" # to perturb all waters in the system
    # waters:
        - "W:15" # chain ID and residue number
        - "W:21"

Metrics
---------

Metrics to track along the simulation

- **atom_dist**: Calculate distance between two atomnumbers. To calculate more than one append them in column as the example below. Default=None

    - The atomdist can be specified via chain:resnum:atomname i.e. A:2:CA

- **rmsd_pdb**: Calculate rmsd of the ligand to a native pdb structure


..  code-block:: yaml

    atom_dist:
        # Distance between the A:2:CA and B:3:CG also between A:5:N and B:3:CG. Append more if desired.
        - "A:2:CA"
        - "B:3:CG"
        - "A:5:N"
        - "B:3:CG"
    rmsd_pdb: "/home/dsoler/native.pdb"


Analysis
-----------

Run a post simulation analysis to extract plots, top poses and clusters.

- **only_analysis**: Analyse PELE simulation without running it.

- **analysis_nclust**: Numbers of clusters out of the simulation, if using the standard clustering method. Default: 10

- **be_column**: Column of the binding energy in the reports starting by 1. Default: 5

- **te_column**: Column of the total energy in the reports starting by 1. Default: 4

- **limit_column**: Specify the column where your external metrics start. Default: 6

- **mae**: To extract the best energy and cluster poses as .mae files with the metrics as properties (schrodinger need it). Default: false

- **analysis**: Whether to run or not the analysis at the end of the simulation. Default: true

- **clustering_method**: If you want to override the default clustering method (Gaussian mixture model), you can set this flag to ``MeanShift`` or ``HDBSCAN``.

- **bandwidth**: Value for the Mean Shift bandwidth (when using the Mean Shift algorithm) or epsilon (when using the HDBSCAN clustering); default = 5.0

- **max_top_clusters**: Maximum number of clusters to be selected. Default = 8.

- **top_clusters_criterion**: Method of selecting top clusters, you can choose one of:

        * "total_25_percentile" - total energy 25th percentile
        * "total_5_percentile" - total energy 5th percentile
        * "total_mean" - total energy mean
        * "total_min" - total energy min
        * "interaction_25_percentile" - interaction energy 25th percentile (default)
        * "interaction_5_percentile" - interaction energy 5th percentile
        * "interaction_mean" - interaction energy mean
        * "interaction_min" - interaction energy min
        * "population" - cluster population

- **cluster_representatives_criterion**: Method of selecting representative structures for each cluster, you can choose one of:

        * "total_25_percentile" - total energy 25th percentile
        * "total_5_percentile" - total energy 5th percentile
        * "total_mean" - total energy mean
        * "total_min" - total energy min
        * "interaction_25_percentile" - interaction energy 25th percentile
        * "interaction_5_percentile" - interaction energy 5th percentile (default)
        * "interaction_mean" - interaction energy mean
        * "interaction_min" - interaction energy min

- **max_top_poses**: Maximum number of top poses to be retrieved. Default = 100.

..  code-block:: yaml

    only_analysis: true
    be_column: 5
    te_column: 4
    limit_column: 6
    mae: true
    clustering_method: "meanshift"
    bandwidth: 7.0
    top_clusters_criterion: "population"
    max_top_clusters: 5
    cluster_representatives_criterion: "interaction_mean"

The bandwidth parameter hugely influences the clustering results, therefore, it might be worth trying out different values depending on your system.
In case of the mean shift algorithm, the bandwidth refers to the maximum RMSD allowed within the cluster, whereas in HDBSCAN to distances between your data points.

Output
----------

Configure the output

- **working_folder**: Name of the main working folder where to store the processed input, control files and the simulation folder. Default="resname_Pele_X" where X is a number.

- **output**: Output folder of the simulation. Default=output

..  code-block:: yaml

    working_folder: "NOR_solvent_OBC"
    output: "output_sim"
