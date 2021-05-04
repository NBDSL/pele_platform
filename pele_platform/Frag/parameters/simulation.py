"""
This module handles the assignment of simulation parameters for FragPELE.
"""
from pele_platform.constants import constants as cs


class FragSimulationParameters(object):
    """
    Base class to assign the simulation parameters for FragPELE.
    """

    def __init__(self, parameters, args):
        """
        Given a Parameters object, it initializes the simulation parameters for
        FragPELE.

        .. todo ::
           * We need to unify all classes that prepare the parameters
             for Frag. We need an abstract class to ensure that all of them
             modify correctly the Parameters object

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        """
        import os

        from pele_platform.constants import constants

        parameters.system = args.frag_core

        # Set number of growing steps
        if args.growing_steps:
            parameters.gr_steps = args.growing_steps
        else:
            parameters.gr_steps = \
                parameters.simulation_params.get("growing_steps", 6)

        # Set number of FragPELE steps
        if args.frag_steps:
            parameters.frag_steps = args.frag_steps
        else:
            parameters.frag_steps = \
                parameters.simulation_params.get("steps_in_gs", 3)

        # Set number of equilibration steps
        if args.frag_eq_steps:
            parameters.frag_eq_steps = args.frag_eq_steps
        else:
            parameters.frag_eq_steps = \
                parameters.simulation_params.get("sampling_steps", 20)

        # Set path to control file
        parameters.control_file = os.path.join(constants.DIR,
                                               "Templates/pele_template.conf")

        # Set FragPELE protocol
        if args.protocol:
            parameters.protocol = args.protocol
        else:
            parameters.simulation_params.get("protocol", "")

        # Set path to topology file
        if parameters.pdb:
            parameters.topology = None
        else:
            parameters.topology = os.path.join("output_pele", "topology.pdb")

        # Set FragPELE restart option
        if args.frag_restart:
            parameters.frag_restart = "-rst"
        else:
            parameters.frag_restart = ""

        # Set trajectory name
        parameters.frag_traj_name = "trajectory"

        # Set use srun parameter for FragPELE
        if parameters.usesrun == "true":
            parameters.usesrun = True
        else:
            parameters.usesrun = False

        parameters.spython = cs.SCHRODINGER
