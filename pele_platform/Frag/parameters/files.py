"""
This module handles the assignment of file parameters for FragPELE.
"""


class FragInputFiles(object):
    """
    Base class to assign the file parameters for FragPELE.
    """

    def __init__(self, parameters, args):
        """
        Given a Parameters object, it initializes the file parameters for
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
        self.input_parameters(parameters, args)
        self.frag_from_core_parameters(parameters, args)
        self.frag_from_sdf_parameters(parameters, args)

    @staticmethod
    def input_parameters(parameters, args):
        """
        Given a Parameters object, it assigns the general parameters
        belonging to the ligand-protein input.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        """
        parameters.input = args.frag_input
        parameters.frag_core_atom = args.frag_core_atom

    @staticmethod
    def frag_from_core_parameters(parameters, args):
        """
        Given a Parameters object, it assigns the parameters that are required
        to run a FragPELE simulation with a series file.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        """
        import os
        import PPP.main as ppp

        parameters.core = args.frag_core

        if args.skip_prep is not None:
            parameters.skip_prep = args.skip_prep
        else:
            parameters.skip_prep = parameters.simulation_params.get("skip_prep", False)

        if not parameters.skip_prep:
            parameters.core_process = \
                os.path.basename(
                    ppp.main(parameters.core, ".", output_pdb=["", ],
                             charge_terminals=args.charge_ter,
                             no_gaps_ter=args.gaps_ter)[0])
        else:
            print("Skipping preprocessing")
            parameters.core_process = parameters.core

        parameters.core_format = parameters.core.rsplit(".")[-1]

    @staticmethod
    def frag_from_sdf_parameters(parameters, args):
        """
        Given a Parameters object, it assigns the parameters that are required
        to run a FragPELE simulation with a sd file.

        Parameters
        ----------
        parameters : a Parameters object
            The Parameters object containing the parameters for PELE
        args : a YamlParser object
            The YamlParser object containing the input parameters chosen
            by the user
        """
        parameters.ligands = args.frag_ligands
        parameters.frag_library = args.frag_library
        parameters.fragment_atom = args.fragment_atom

