"""
This module holds the classes and methods to build and manage the parameters
that are required by the PELE Platform to run the different workflows.
"""

__all__ = ["ParametersBuilder", "Parameters"]


from pele_platform.Utilities.Parameters.SimulationParams import \
    simulation_params
from pele_platform.Utilities.Parameters.SimulationFolders import \
    simulation_folders


class ParametersBuilder(object):
    """
    It builds a Parameters instance and it creates the required folders and
    files
    """

    def __init__(self):
        """
        It initializes a ParametersBuilder object.
        """
        self._parameters = None
        self._initialized = False

    def build_adaptive_variables(self, args):
        """
        It builds the parameters for adaptive, according to the arguments
        that are supplied, and returns the corresponding Parameters
        instance.

        Parameters
        ----------
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user

        Returns
        -------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        """
        import os
        from pele_platform.constants import constants
        from pele_platform.features import adaptive
        from pele_platform.Utilities.Helpers import helpers

        # Define main PELE directory
        main_dir = os.path.abspath("{}_Pele".format(args.residue))

        # Set the PELE directory
        # In case that folder is not set by the user, we will try to suggest
        # the best candidate considering whether we want to restart a previous
        # simulation or we want to run a new one from scratch.
        if not args.folder:
            # Check if the simulation is being restarted or not
            # TODO the restart flag is undocumented and counterintuitive
            #      since it apparently is doing the opposite of
            #      adaptive_restart
            if args.restart in constants.FIRST_RESTART:
                pele_dir = helpers.get_next_peledir(main_dir)
            else:
                pele_dir = helpers.get_latest_peledir(main_dir)

            # Check if the adaptive simulation is being restarted or not
            # Also check if we are only running the analysis
            if args.adaptive_restart or args.only_analysis:
                # Take the last folder name
                pele_dir = helpers.get_latest_peledir(main_dir)
            else:
                # Get a new folder name
                pele_dir = helpers.get_next_peledir(main_dir)

        # In case that the user has specified the output folder, we will
        # use it, regardless it already exists.
        else:
            pele_dir = os.path.abspath(args.folder)

        # Retrieve the specific args for adaptive
        specific_args = adaptive.retrieve_software_settings(args,
                                                            pele_dir)

        # Add pele_dir
        specific_args['pele_dir'] = pele_dir

        # Initialize Parameters object
        self._parameters = Parameters(args, specific_args)
        self._initialized = True

        # Set software
        self.parameters.software = "Adaptive"

        return self.parameters

    def build_frag_variables(self, args):
        """
        It builds the parameters for frag, according to the arguments
        that are supplied, and returns the corresponding Parameters
        instance.

        Parameters
        ----------
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user

        Returns
        -------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        """
        from pele_platform.features import frag
        from pele_platform.Frag.parameters import FragWaterParams
        from pele_platform.Frag.parameters import FragInputFiles
        from pele_platform.Frag.parameters import FragSimulationParameters
        from pele_platform.Frag.parameters import FragMetrics
        from pele_platform.Frag.parameters import FragOptionalParameters

        # Retrieve the specific args for FragPELE
        specific_args = frag.retrieve_software_settings(args)

        # Initialize Parameters object
        self._parameters = Parameters(args, specific_args,
                                      initialize_simulation_paths=False)
        self._initialized = True

        # Initialize water parameters for FragPELE
        FragWaterParams(self.parameters, args)

        # Initialize file parameters for FragPELE
        FragInputFiles(self.parameters, args)

        # Initialize simulation parameters for FragPELE
        FragSimulationParameters(self.parameters, args)

        # Initialize metrics parameters for FragPELE
        FragMetrics(self.parameters, args)

        # Initialize optional parameters for FragPELE
        FragOptionalParameters(self.parameters, args)

        # Keep initial arguments
        self.parameters.args = args

        # Create logger
        self.parameters.create_logger(".")

        # Set software
        self.parameters.software = "Frag"

        return self.parameters

    @property
    def initialized(self):
        """
        It returns the initialization state of the builder. If this builder
        has already built a Parameters object before, it will return a True.
        Otherwise, it will return a False.

        Returns
        -------
        initialized : bool
            Whether this ParametersBuilder has already been used to build
            a Parameters object or not
        """
        return self._initialized

    @property
    def parameters(self):
        """
        It returns the Parameters object that has been built by the current
        ParametersBuilder instance, if any. In case that it has not been built
        yet (ParametersBuilder not initialized), it will return None.

        Returns
        -------
        _parameters : a Parameters object
            The Parameters object containing the right parameters according
            to the input given by the user
        """
        return self._parameters


class Parameters(simulation_params.SimulationParams,
                 simulation_folders.SimulationPaths):
    """
    Base parameters class where the PELE Platform parameters are stored
    and manipulated.

    .. todo ::
       * This class should be serializable. For example, we should be able to
         represent it as a json format at any time. We need to implement
         an abstract class for this purpose, something like in:
         https://github.com/openforcefield/openff-toolkit/blob/master/openff/toolkit/utils/serialization.py#L25

        * The args parameter that is required must be refactored. Among other
          things, it must support a Python API and its initialization and
          modification must be as easy and straightforward as possible.

       * Consider moving the functionality about the creation of files
         and folders outside this method.
    """

    def __init__(self, args, specific_args=None,
                 initialize_simulation_params=True,
                 initialize_simulation_paths=True):
        """
        It initializes a Parameters object.

        .. todo ::
           * We need to refactor the inheritance from
              - simulation_params.SimulationParams and
              - simulation_folders.SimulationPaths
             since the parameters initialization is currently very
             complex and confusing.

        Parameters
        ----------
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        specific_args : dict
            The dictionary containing the specific parameters for this
            Parameter object
        initialize_simulation_params : bool
            Whether to initialize the simulation parameters or not. Default
            is True
        initialize_simulation_paths : bool
            Whether to initialize the simulation paths or not. Default
            is True
        """
        # CA interval is first initialized to None
        self.ca_interval = None

        # If specific_args is not defined, assign and empty dictionary
        if specific_args is None:
            specific_args = dict()

        # Set specific parameters, they need to be set before initializing
        # the rest
        for key, value in specific_args.items():
            setattr(self, key, value)

        # Initialize the parameters from parent classes
        if initialize_simulation_params:
            simulation_params.SimulationParams.__init__(self, args)
        if initialize_simulation_paths:
            simulation_folders.SimulationPaths.__init__(self, args)

        # We need to set the specific arguments again
        for key, value in specific_args.items():
            setattr(self, key, value)

    def create_files_and_folders(self):
        if not self.adaptive_restart and not self.only_analysis:
            self.create_folders()
            self.create_files()
            self.create_logger()
        else:
            self.create_logger()

    def create_folders(self):
        """
            Create pele folders
        """
        from pele_platform.Utilities.Helpers import helpers

        for folder in self.folders:
            helpers.create_dir(self.pele_dir, folder)

    def create_files(self):
        """
            Copy templates
        """
        import os
        import shutil

        # Actions
        for file, destination_name in zip(self.files, self.file_names):
            shutil.copy(file, os.path.join(self.pele_dir, destination_name))

    def create_logger(self, directory=None):
        """
        Given a directory, it initiates a logger in that place.

        Parameters
        ----------
        directory : str
            The directory where the logger will be stored. Default uses the
            main PELE directory
        """
        import os
        import logging

        directory = directory if directory else self.pele_dir
        log_name = os.path.join(directory, "{}.log".format(self.residue))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
        if not self.adaptive_restart:
            file_handler = logging.FileHandler(log_name, mode='w')
        else:
            file_handler = logging.FileHandler(log_name, mode='a')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def to_json(self, indent=None):
        """
        Return a JSON serialized representation.
        Specification: https://www.json.org/

    .. todo ::
       * Consider moving this method to a general abstract class to
         serialize classes.

        Parameters
        ----------
        indent : int, optional, default=None
            If not None, will pretty-print with specified number of spaces
            for indentation

        Returns
        -------
        serialized : str
            A JSON serialized representation of the object
        """
        import json

        d = self.to_dict()

        return json.dumps(d, indent=indent)

    def to_dict(self):
        """
        Return a dictionary representation of the Parameters object.

        Returns
        -------
        parameters_dict : an OrderedDict object
            The dictionary representation of the current Parameters object.
        """
        from collections import OrderedDict

        parameters_dict = OrderedDict()

        for attribute, value in self.__dict__.items():
            # Skip special and private attributes
            if attribute.startswith('__') and attribute.endswith('__'):
                continue

            parameters_dict[attribute] = value

        return parameters_dict
