import os
import shutil
import warnings
import random
import logging
import MSM_PELE.constants as cs


def set_pele_env(folders, files, forcefield, pele_dir):
    pele_env = Pele_env_Builder(folders, files, forcefield, pele_dir)
    pele_env.folder_levels()
    pele_env.file_dist()


class Pele_env_Builder(object):
    """
        Base class wher the needed pele environment
        is build by creating folders and files
    """

    def __init__(self, folders, files, forcefield, system, residue, cpus, restart, native, chain, mae_lig):
        self.folders = folders
        self.files = files
        self.system = system
        self.forcefield = forcefield
        self.residue = residue
        self.templates = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), "PeleTemplates"))
        self.cpus = cpus
        self.restart = restart
        self.native = native
        self.chain = chain
        self.mae_lig = mae_lig

        self.build_constant_paths()

    def build_constant_paths(self):

        self.template = None
        self.rotamers_file = None
        self.random_num = random.randrange(1, 70000)
        self.license = '''"{}"'''.format(cs.LICENSE)
        self.equil_steps = int(cs.EQ_STEPS/self.cpus) if self.cpus < cs.EQ_STEPS else 1
        pele_dir = os.path.abspath("{}_Pele".format(self.residue))
        self.pele_dir = is_repited(pele_dir) if self.restart == "all" else is_last(pele_dir)
        self.native = cs.NATIVE.format(os.path.abspath(self.native), self.chain) if self.native else self.native
        if self.mae_lig:
            self.system_fix = os.path.join(self.pele_dir, "{}_complex_processed.pdb".format(os.path.abspath(os.path.splitext(self.system)[0])))
        else:
            self.system_fix = os.path.join(self.pele_dir, "{}_processed.pdb".format(os.path.abspath(os.path.splitext(self.system)[0])))
        self.adap_ex_input = os.path.join(self.pele_dir, os.path.basename(self.system_fix))
        self.adap_ex_output = os.path.join(self.pele_dir, "output_adaptive_exit")
        self.cluster_output = os.path.join(self.pele_dir, "output_clustering")
        self.adap_l_input = "{}/initial_*"
        self.adap_l_output = os.path.join(self.pele_dir, "output_pele")
        self.ad_ex_temp = os.path.join(self.pele_dir, "adaptive_exit.conf")
        self.ad_l_temp = os.path.join(self.pele_dir, "adaptive_long.conf")
        self.pele_exit_temp = os.path.join(self.pele_dir, "pele_exit.conf")
        self.pele_temp = os.path.join(self.pele_dir, "pele.conf")
        self.box_temp = os.path.join(self.pele_dir, "box.pdb")
        self.clusters_output = os.path.join(self.cluster_output, "clusters_40_KMeans_allSnapshots.pdb")
        self.ligand_ref = os.path.join(self.pele_dir, "ligand.pdb")
       

    def create(self):
        if self.restart == "all":
			self.create_folders()
			self.create_files()
			self.create_logger()
        else:
            self.create_logger()

    def create_folders(self):
        """
            Create pele folders
        """

        for folder in self.folders:
            self.create_dir(self.pele_dir, folder)

    def create_files(self):
        """
            Copy templates
        """

        # Actions
        for file in self.files:
            self.copy(file, self.pele_dir)


    def create_dir(self, base_dir, extension=None):
        """
            Class Method to manage
            directory creation only if that
            ones doesn't exist

            Location:
                base_dir+extension
                or base_dir if extension is None
        """
        if extension:
            path = os.path.join(base_dir, extension)
            if os.path.isdir(path):
                warnings.warn("Directory {} already exists.".format(path), RuntimeWarning)
            else:
                os.makedirs(path)
        else:
            if os.path.isdir(base_dir):
                warnings.warn("Directory {} already exists.".format(base_dir), RuntimeWarning)
            else:
                os.makedirs(base_dir)

    def copy(self, standard, destination, user=None):
        if user:
            shutil.copy(user, os.path.join(self.pele_dir, standard))
        else:
            shutil.copy(standard, self.pele_dir)
        return os.path.join(self.pele_dir, standard)

    def create_logger(self):
        log_name = os.path.join(self.pele_dir, "{}.log".format(self.residue))
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s:%(levelname)s:%(message)s")
        if self.restart == "all":
			file_handler = logging.FileHandler(log_name, mode='w')
        else:
			file_handler = logging.FileHandler(log_name, mode='a')
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)



def is_repited(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
        if chunk != "Pele":
            if original_dir:
                original_dir = "{}_{}".format(original_dir, chunk)
            else:
                original_dir = chunk
        else:
            break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
        i = 1
    if os.path.isdir(pele_dir):
		new_pele_dir = "{}_Pele_{}".format(original_dir, i)
		new_pele_dir = is_repited(new_pele_dir)
		return new_pele_dir
    else:
		return pele_dir

def is_last(pele_dir):

    original_dir = None
    split_dir = pele_dir.split("_")
    for chunk in split_dir:
		if chunk != "Pele":
			if original_dir:
 				original_dir = "{}_{}".format(original_dir, chunk)
			else:
				original_dir = chunk
		else:
			break
    if split_dir[-1].isdigit():
        i = split_dir[-1]
        i = int(i) + 1 
    else:
		i = 1 

    if os.path.isdir(pele_dir):
            new_pele_dir = "{}_Pele_{}".format(original_dir, i)
            if not os.path.isdir(new_pele_dir):
                return pele_dir
            else:
			    new_pele_dir = is_last(new_pele_dir)
			    return new_pele_dir
    else:
        return pele_dir

