import glob
import os
import subprocess
import shutil

from rdkit import Chem
from pele_platform.constants import constants as cs

OUTPUT = "input.conf"


def get_symmetry_groups(mol):
    """
    Computes the symmetry class for each atom and returns a list with the idx of non-symmetric atoms.
    """
    rank = {} 
    symmetry_list = []
    symmetry_rank_list = []
    counter = 0
    
    for atom in mol.GetAtoms():
        rank[atom.GetIdx()] = list(Chem.CanonicalRankAtoms(mol, breakTies=False))[counter]
        counter += 1
    
    for idx, symmetry_rank in rank.items():
        if symmetry_rank not in symmetry_rank_list:
            symmetry_rank_list.append(symmetry_rank)
            symmetry_list.append(idx)
    return symmetry_list


def growing_sites(fragment,
                  user_bond):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.
    Output - list of strings representing sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    bonds = []
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    symmetry_list = get_symmetry_groups(mol)
    if mol:
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H" and n.GetIdx() in symmetry_list]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                h_name = h.GetMonomerInfo().GetName().strip()
                bonds.append("{} {} {}-{}".format(fragment, user_bond, at_name, h_name))
    return bonds


def sdf_to_pdb(file_list,
               path, logger, tmpdirname):

    out = []
    if file_list:
        converted_mae = []
        output = []

        # convert all SDF to MAE
        schrodinger_path = os.path.join(cs.SCHRODINGER, "utilities/structconvert")
        command_mae = "{} -isd {} -omae {}"
        command_pdb = "{} -imae {} -opdb {}"
        for f in file_list:
            shutil.copy(f, tmpdirname)
            fout = os.path.splitext(os.path.basename(f))[0] + ".mae"
            fout_path = '%s/%s' % (tmpdirname, os.path.basename(f))
            try:
                command_mae = command_mae.format(schrodinger_path, fout_path, fout)
                subprocess.call(command_mae.split())
                converted_mae.append(fout)
            except Exception as e:
                logger.info("Error occured while converting SD files to mae.", e)
           
        # convert all MAE to PDB, it will result in a lot of numbered pdb files
        for c in converted_mae:
            shutil.move(c, tmpdirname)
            c = tmpdirname + '/' + c
            fout = c.replace(".mae", ".pdb")
            
            try:
                command_pdb = command_pdb.format(schrodinger_path, c, fout)
                subprocess.call(command_pdb.split())
                os.remove(c)
            except Exception as e:
                logger.info("Error occured while converting mae to PDB.", e)
        
        pdb_pattern = '%s/%s' % (tmpdirname, converted_mae[0])
        converted_pdb = glob.glob(pdb_pattern[:-4]+"*"+".pdb")
        # ~~~ If it's stupid but it works (?), it isn't stupid. ~~~
        
        # read in PDB file created by Schrodinger, substitute residue name and add chain ID
        for c in converted_pdb:
            with open(c, "r") as fin:
                lines = fin.readlines()
                new_lines = []
                for line in lines:
                    if line.startswith("HETATM") or line.startswith("ATOM"):
                        new_lines.append(line)
            
            new_lines = [l.replace("UNK", "GRW") for l in new_lines if "UNK" in l]
            new_lines = [l[:21]+"L"+l[22:] for l in new_lines]

            with open(c, "r+") as fout:
                for line in new_lines:
                    fout.write(line)
        out = converted_pdb
    return out


def get_library(frag_library):
    directory = os.path.dirname(os.path.abspath(__file__))                                                                                                                              
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip())                                                                 
    if not os.path.exists(path):                                                                                                                                                        
        raise OSError(f"File {frag_library} doesn't exist and is not one of our internal libraries. Please check the frag_library flag in input.yaml.")   
    return path


def get_fragment_files(path, logger, tmpdirname):

    fragment_files = []                                                                                                                                                                 
    extensions = ['*.pdb', '*.sdf']                                                                                                                                                     
    
    for e in extensions:                                                                                                                                                                
        fragment_files.extend(glob.glob(os.path.join(path, e.upper())))                                                                                                                 
        fragment_files.extend(glob.glob(os.path.join(path, e.lower())))                                                                                                                 
                                                                                                                                                                                        
    # convert SDF to PDB, if necessary                                                                                                                                                  
    sdf_files = [elem for elem in fragment_files if ".sdf" in elem.lower()]                                                                                                             
    pdb_files = [elem for elem in fragment_files if ".pdb" in elem.lower()]                                                                                                             
    all_files = pdb_files + sdf_to_pdb(sdf_files, path, logger, tmpdirname)                                                                                                                                                                               
    return all_files


def write_config_file(output_name,
                      bond_list):

    with open(output_name, "w+") as conf_file:
        for line in bond_list:
            conf_file.write(line+"\n")


def main(user_bond, frag_library, logger, tmpdirname):

    # find the library and extract fragments
    path = get_library(frag_library)
    all_files = get_fragment_files(path, logger, tmpdirname) 
    
    # get all possible growing sites
    bond_list = []
    for file in all_files:
        bond_list.extend(growing_sites(file, user_bond))
    
    # write input.conf 
    write_config_file(OUTPUT, bond_list)
    
    return OUTPUT
