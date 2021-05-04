import pandas as pd
import glob
import os
import shutil
from pele_platform.PPI.preparation import prepare_structure
from pele_platform.constants import constants as cs
from pele_platform import main


test_path = os.path.join(cs.DIR, "Examples")



yaml = os.path.join(test_path, "PPI/input_skipref.yaml")
def test_ppi_skipref(yaml=yaml):

    #Function to test
    job, _ = main.run_platform_from_yaml(yaml)

    # checkpoints
    files_refinement = glob.glob(os.path.join(job.pele_dir, "refinement_simulation/results/top_poses/*.pdb"))

    # test
    assert not files_refinement

yaml = os.path.join(test_path, "PPI/input.yaml")
def test_ppi_default(energy_result=-1.99, yaml=yaml):
  
    #Function to test
    job, job2 = main.run_platform_from_yaml(yaml)

    # checkpoints
    output_csv = pd.read_csv(os.path.join(job.pele_dir, "output/clustering_output.csv"))
    best_energy = round(output_csv["binding_energy"].min(), 2)
    nfiles = len(glob.glob(os.path.join(os.path.dirname(job.pele_dir), "refinement_input/*.pdb")))
    nfiles_refinement = len(glob.glob(os.path.join(job2.pele_dir, "results/top_poses/*.pdb")))

    # test
    assert nfiles == job.n_components 
    assert best_energy == energy_result
    assert nfiles_refinement

def test_prepare_structure():

    protein_file = os.path.join(test_path, "PPI/1tnf_prep.pdb")
    new_protein_file = "1tnf_prep_prep.pdb"
    ligand_pdb = os.path.join(test_path, "PPI/1tnf_ligand.pdb")
    chain = ["A", "B"]

    prepare_structure(protein_file, ligand_pdb, chain, remove_water=True)

    with open(new_protein_file, "r") as file:
        lines = file.readlines()
        chains = []
        no_water = True
        for line in lines:
            if "HOH" in line:
                no_water = False
            if line.startswith("HETATM") or line.startswith("ATOM"):
                chains.append(line[21:22].strip())
    
    assert "C" not in chains
    assert no_water
    
    for f in glob.glob("*_prep.pdb"):
        os.remove(f)


yaml = os.path.join(test_path, "PPI/input_skipref.yaml")

def test_ppi_skipref(energy_result=-2.18, yaml=yaml):

    #Function to test
    job, _ = main.run_platform_from_yaml(yaml)

    # checkpoints
    files_refinement = glob.glob(os.path.join(job.pele_dir, "refinement_simulation/results/top_poses/*.pdb"))

    # test
    assert not files_refinement


yaml = os.path.join(test_path, "PPI/input_folder.yaml")
def test_working_folder(yaml=yaml, output="ppi_folder"):
    if os.path.exists(output): shutil.rmtree(output)
    job, _ = main.run_platform_from_yaml(yaml)
    assert os.path.exists(job.folder)
