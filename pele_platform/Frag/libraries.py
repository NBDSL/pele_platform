import glob
import os
import subprocess
import shutil
import sys

from rdkit import Chem
from pele_platform.constants import constants as cs

OUTPUT = "input.conf"

def get_symmetry_groups(mol):
    """
    Computes the symmetry class for each atom and returns a list with the idx of non-symmetric atoms.

    Parameters
    ----------
    mol : rdkit molecule object.
          Fragment from custom-made library.

    Returns
    -------
    symmetry_list : list
                    List with atom indices.
    """
    rank = {} 
    symmetry_list = []
    symmetry_rank_list = []
    counter = 0
    
    for counter, atom in enumerate(mol.GetAtoms()):
        rank[atom.GetIdx()] = list(Chem.CanonicalRankAtoms(mol, breakTies=False))[counter]
    for idx, symmetry_rank in rank.items():
        if symmetry_rank not in symmetry_rank_list:
            symmetry_rank_list.append(symmetry_rank)
            symmetry_list.append(idx)
    return symmetry_list


def growing_sites(fragment,
                  user_bond,
                  fragment_atom):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.

    Parameters
    ----------
    fragment : string
               Path to fragment pdb file.
    
    user_bond : string
                Connection point from which the user wants to grow the fragments.

    Returns
    -------
    bonds : list
            List of strings representing sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    bonds = []
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    symmetry_list = get_symmetry_groups(mol)
    if fragment_atom is not None and mol:
        # THE FRAGMENT_ATOM SHOULD BE AN HYDROGEN FROM WHERE WE WANT TO ATTACH THE FRAGMENT
        for atom in mol.GetAtoms():
            if atom.GetPDBResidueInfo().GetName().strip() == fragment_atom.strip():
                fragment_atom_idx = atom.GetIdx()
        fragment_bonds = [(x.GetBeginAtomIdx(), x.GetEndAtomIdx()) for x in mol.GetBonds()]
        for bond in fragment_bonds:
            if fragment_atom_idx in bond:
                if bond[0] == fragment_atom_idx:
                    core_id = bond[1]
                else:
                    core_id = bond[0]
        for a in mol.GetAtoms():
            if a.GetIdx() == core_id:
                core_id_name = a.GetPDBResidueInfo().GetName().strip()
        bonds.append("{} {} {}-{}".format(fragment, user_bond, core_id_name, fragment_atom.strip()))
    else:
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
               logger, tmpdirname):
    """
    Converts sdf files to pdb.

    Parameters
    ----------
    file_list : list
                List of paths of fragments in sdf format.
    
    logger : string
             File with status messages
    
    tmpdirname : string
                 Path of temporary directory.

    Returns
    -------
    out : list
          List of paths with converted sdf files.
    """

    out = []
    if file_list:
        converted_mae = []

        # convert all SDF to MAE
        schrodinger_path = os.path.join(cs.SCHRODINGER, "utilities/structconvert")
        command_mae = "{} -isd {} -omae {}"
        command_pdb = "{} -imae {} -opdb {}"
        for file in file_list:
            shutil.copy(file, tmpdirname)
            fout = os.path.splitext(os.path.basename(file))[0] + ".mae"
            fout_path = os.path.join(tmpdirname, os.path.basename(file))
            try:
                command_mae = command_mae.format(schrodinger_path, fout_path, fout)
                subprocess.call(command_mae.split())
                converted_mae.append(fout)
            except Exception as e:
                logger.info("Error occured while converting SD files to mae.", e)
        # convert all MAE to PDB, it will result in a lot of numbered pdb files
        for c in converted_mae:
            shutil.move(c, tmpdirname)
            c = os.path.join(tmpdirname, c)
            fout = c.replace(".mae", ".pdb")

            try:
                command_pdb = command_pdb.format(schrodinger_path, c, fout)
                subprocess.call(command_pdb.split())
                os.remove(c)
            except Exception as e:
                logger.info("Error occured while converting mae to PDB.", e)
        

        pdb_pattern = os.path.join(tmpdirname, converted_mae[0])
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

            with open(c, "w") as fout:
                for line in new_lines:
                    fout.write(line)
        out = converted_pdb
    return out


def get_library(frag_library):
    """
    Checks the path of the fragment library provided on the input.yaml file.

    Parameters
    ----------
    frag_library : string
                   Path to fragment library.

    Returns
    -------
    path : string
           Path to the fragment library.
    """
    directory = os.path.dirname(os.path.abspath(__file__))                                                                                                                              
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip())                                                                 
    if not os.path.exists(path):                                                                                                                                                        
        raise OSError(f"File {frag_library} doesn't exist and is not one of our internal libraries. Please check the frag_library flag in input.yaml.")   
    return path


def get_fragment_files(path,
                       logger, tmpdirname):
    """
    Gets all pdb and sdf files of each fragment in the library.

    Parameters
    ----------
    path : string
           Path to the fragment library.

    logger : string
             File with status messages.

    tmpdirname : string
                 Path of temporary directory.

    Returns
    -------
    all_files : list
                List of paths of the fragments in the fragment library.
    """
    fragment_files = []                                                                                                                                                                 
    extensions = ['*.pdb', '*.sdf']                                                                                                                                                     
    
    for e in extensions:                                                                                                                                                                
        fragment_files.extend(glob.glob(os.path.join(path, e.upper())))                                                                                                                 
        fragment_files.extend(glob.glob(os.path.join(path, e.lower())))                                                                                                                 
                                                                                                                                                                                        
    # convert SDF to PDB, if necessary                                                                                                                                                  
    sdf_files = [elem for elem in fragment_files if ".sdf" in elem.lower()]                                                                                                             
    pdb_files = [elem for elem in fragment_files if ".pdb" in elem.lower()]
    all_files = pdb_files + sdf_to_pdb(sdf_files, logger, tmpdirname)
    return all_files


def write_config_file(output_name,
                      bond_list, frag_restart, frag_core):
    """
    Generates the configuration file.
    """
    if frag_restart:
        for i in bond_list:

            name= os.getcwd() + "/"+ frag_core.split("/")[-1][:-4]+"_processed_"+ \
                  i.split("/")[6].split(" ")[0].replace(".pdb","") + \
                  i.split("/")[6].split(" ")[1].replace(" ","") + \
                  i.split("/")[6].split(" ")[2].replace(" ","")  + \
                  "/" + i.split("/")[6].split(" ")[0].replace(".pdb","") + \
                  i.split("/")[6].split(" ")[1].replace(" ", "") + \
                  i.split("/")[6].split(" ")[2].replace(" ", "") + "_top.pdb"

            if os.path.isfile(name):
                bond_list.remove(i)
    with open(output_name, "w+") as conf_file:
        for line in bond_list:
            conf_file.write(line+"\n")


def main(frag_core,
         user_bond,
         frag_library,
         logger,
         fragment_atom,
         frag_restart,
         tmpdirname):
    # find the library and extract fragments
    path = get_library(frag_library)
    all_files = get_fragment_files(path, logger, tmpdirname) 
    
    # get all possible growing sites
    bond_list = []
    for file in all_files:
        bond_list.extend(growing_sites(file, user_bond, fragment_atom))
    
    # write input.conf 
    write_config_file(OUTPUT, bond_list, frag_restart, frag_core)
    
    return OUTPUT
