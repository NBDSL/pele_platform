# -*- coding: utf-8 -*-
"""
This is the main module and it is designed to run the PELE platform from
the command-line.
"""

__author__ = "Nostrum Biodiscovery"
__email__ = "pelesupport@nostrumbiodiscovery.com"
__license__ = "Apache-2.0"


def parse_args(args=[]):
    """
    Command line parser.

    Parameters
    ----------
    args : list[str]
        The list of strings to parse. Default is an empty list. If empty,
        input arguments are retrieved from the command-line

    Returns
    -------
    input_yaml : str
        The path pointing to the input yaml file
    """
    # Parser setup
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Automatic platform to launch PELE simulations')
    parser.add_argument('input_file', type=str, help='Yaml input file')

    # Parse arguments
    parsed_args = parser.parse_args(args) if args else parser.parse_args()

    # Extract yaml file (the only argument we need to retrieve)
    input_yaml = parsed_args.input_file

    return input_yaml


def run_platform_from_yaml(input_yaml):
    """
    High level function to run the PELE platform from an input yaml file.

    It will:
    1) Parse the input.yaml
    2) Launch job
    3) Return job parameters

    Parameters
    ----------
    input_yaml : str
        The path pointing to the input yaml file

    Returns
    -------
    job_params : a Parameters object
        The corresponding Parameters object with the parameters of the
        simulation
    """
    # Generate the yaml object from the input yaml
    from pele_platform.Utilities.Helpers import yaml_parser
    from pele_platform.Checker import valid_flags

    yaml_obj = yaml_parser.YamlParser(input_yaml,
                                      valid_flags.VALID_FLAGS_PLATFORM)

    # Attempt to parse the yaml object
    try:
        yaml_obj.read()
    except AttributeError:
        from pele_platform.Errors import custom_errors

        raise custom_errors.WrongYamlFile(
            "Input file: {}".format(input_yaml)
            + " does not look like a correct yaml file")

    # Initialize job launcher
    from pele_platform.Utilities.Helpers.launcher import Launcher

    launcher = Launcher(yaml_obj)

    # Run launcher
    job_params = launcher.launch()

    # Return job parameters
    return job_params


# Main workflow to be executed
if __name__ == "__main__":
    # Avoid backend issues in matplotlib
    import matplotlib

    matplotlib.use("Agg")

    # Avoid Python2
    from .Checker.python_version import check_python_version

    check_python_version()

    # Parse yaml file from command-line arguments
    yaml_file = parse_args()

    # Call platform runner
    job_params = run_platform_from_yaml(yaml_file)
