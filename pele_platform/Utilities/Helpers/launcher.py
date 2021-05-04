from dataclasses import dataclass
import pele_platform.Checker.main as ck
import pele_platform.Frag.simulation as fr
import pele_platform.Adaptive.simulation as ad
import pele_platform.site_finder.main as al
import pele_platform.gpcr.main as gpcr
import pele_platform.out_in.main as outin
from pele_platform.PPI.main import run_ppi
from pele_platform.enzyme_engineering.saturated_mutagenesis import SaturatedMutagenesis
import pele_platform.Utilities.Parameters.parameters as pv
import argparse


@dataclass
class Launcher:

    _args: argparse.ArgumentParser
    frag: str = "frag"
    ppi: str = "PPI"
    site_finder: str = "site_finder"
    gpcr_orth: str = "gpcr_orth"
    out_in: str = "out_in"
    adaptive: str = "adaptive"
    saturated_mutagenesis: str = "saturated_mutagenesis"
    interaction_restrictions: str = "interaction_restrictions"

    def launch(self) -> pv.ParametersBuilder:
        # Launch package from input.yaml
        self._define_package_to_run()
        job_variables = self.launch_package(self._args.package, no_check=self._args.no_check)
        return job_variables

    def launch_package(self, package: str, no_check=False) -> pv.ParametersBuilder:
        # Launch package from API
        if not no_check:
            ck.check_executable_and_env_variables(self._args)
        if package == self.adaptive or package == "interaction_restrictions":
            job_variables = ad.run_adaptive(self._args)
        elif package == self.gpcr_orth:
            job_variables = gpcr.GpcrLauncher(self._args).run_gpcr_simulation()
        elif package == self.out_in:
            job_variables = outin.OutInLauncher(self._args).run_outin_simulation()
        elif package == self.site_finder:
            job_variables = al.SiteFinderLauncher(self._args).run_site_finder()
        elif package == self.ppi:
            job_variables = run_ppi(self._args)
        elif package == self.saturated_mutagenesis:
            job_variables = SaturatedMutagenesis(self._args).run()
        elif package == self.frag:
            # Set variables and input ready
            job_variables = fr.FragRunner(self._args).run_simulation()
        return job_variables

    def _define_package_to_run(self) -> None:
        # Define package being run from input.yaml flags
        if self._args.frag_core:
            self._args.package = self.frag
        elif self._args.ppi:
            self._args.package = self.ppi
        elif self._args.site_finder:
            self._args.package = self.site_finder
        elif self._args.gpcr_orth:
            self._args.package = self.gpcr_orth
        elif self._args.out_in:
            self._args.package = self.out_in
        elif self._args.saturated_mutagenesis:
            self._args.package = self.saturated_mutagenesis
        elif self._args.interaction_restrictions:
            self._args.package = self.interaction_restrictions
        else: 
            self._args.package = self.adaptive
