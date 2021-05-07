from dataclasses import dataclass, field
from typing import List
from Bio.PDB import PDBParser, PDBIO, Selection, NeighborSearch, Vector
#import xpdb  # this is the module described below
import glob
import numpy as np
import os
import pele_platform.constants.constants as cs
import pele_platform.Errors.custom_errors as ce
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.pele_params as pp
import sys

TEMPLATE = '{{"watersToPerturb": {{"links": {{"ids": [{index}] }}}}, "Box": {{"radius": {radius}, "fixedCenter": [{com}], "type": "sphericalBox"}}}}'

@dataclass
class WaterIncluder():

    input_pdbs: list
    n_waters: int
    user_waters: List=field(default_factory=lambda: [])
    ligand_perturbation_params: str=""
    ligand_residue: str=""
    water_center: bool=False
    water_radius: bool=False
    water_to_exclude: List=field(default_factory=lambda: [])
    sim_path: str="."
    allow_empty_selectors: bool=False
    water_temp: int=5000
    water_trials: int=10000
    water_overlap: float=0.78
    water_constr: float=0
    water_freq: int=1
    all_waters: str=""
    perturbed_waters=[]
    test: bool=False

    
    def run(self):
        self.set_empty_selectors()
        self.set_user_waters()
        self.water_checker()
        self.add_water()
        self.set_water_control_file()
    
    def set_empty_selectors(self):
        #If True include in control file else null
        self.allow_empty_selectors = '"allowEmptyWaterSelectors": true,' if self.allow_empty_selectors else ""

    def set_user_waters(self):
        #If all_waters extract them from the any pdb
        if self.user_waters == "all_waters":
            self.user_waters = hp.retrieve_all_waters(self.input_pdbs[0])
        
        
    def set_parameters(self):
        if self.ligand_perturbation_params:
            self.ligand_perturbation_params = self.ligand_perturbation_params.rstrip("]\n") + pp.WATER_PARAMS
        else:
            self.ligand_perturbation_params = ', "parametersChanges" : [ ' + pp.WATER_PARAMS.replace(",{", "{")

    def set_box_center(self, inp=False):
        if self.ligand_residue:
            com = ligand_com(inp, self.ligand_residue)
            com_format = ['{:.10f}'.format(elem) for elem in com[0]]
            self.water_center = ", ".join(com_format)
        else:
            com_format = ['{:.10f}'.format(coord) for coord in self.water_center]
            self.water_center = ", ".join(com_format)


    def set_box_radius(self):
        if self.water_radius:
            self.water_radius =  self.water_radius
        else:
            self.water_radius = 6

    def set_water_input(self, inp):
        try:
            self.set_box_center(inp)
        except TypeError:
            raise ce.NotCenterOfWaterBox('Center of water box not found. Specify with the next flag\n\n"water_center":\n - 18\n - 17\n - 18')
        self.set_box_radius()
        waters = sorted(hp.retrieve_all_waters(inp, exclude=self.water_to_exclude))
        waters = ", ".join(['"{}"'.format(water) for water in waters])
        self.all_waters += waters + ", "
        return TEMPLATE.format(index=waters, radius=self.water_radius, com=self.water_center)
             

    def set_water_control_file(self):
        if self.n_waters != 0 or self.user_waters:
            self.set_parameters()
            water_string = [self.set_water_input(inp).strip("'") for inp in self.input_pdbs]
            self.all_waters = self.all_waters.rstrip(", ")
            self.water_line = cs.WATER.format(self.all_waters, self.allow_empty_selectors, self.water_temp, self.water_trials,
                                        self.water_overlap, self.water_constr, water_string).replace("'", "")
            self.water_energy = None  # TEMPORARY FIX
    
        else:
            self.water_energy = None
            self.water_radius = None
            self.water_center = None
            self.water_line = ""
    
    
    def water_checker(self):
        max_water = 4
        min_water = 1
        if self.n_waters:
            if self.n_waters > max_water:
                raise ValueError("Maximum {} water molecules are allowed.".format(max_water))
            elif self.n_waters < min_water:
                raise ValueError(
                    "Number of water molecules (n_waters) has to be between {} and {}".format(min_water, max_water))
    
    
    def add_water(self):
        if self.test:
            np.random.seed(42)
    
    
        output = []
        n_inputs = len(self.input_pdbs)
        water_coords = []
        resnums = []
        atomnums = []
        chains = []
        resnames = []
    
        # get maximum residue and atom numbers keep original waters
        with open(self.input_pdbs[0], "r") as file:
            lines = file.readlines()
            conect = [line for line in lines if "CONECT" in line]
            pdb_lines = [line for line in lines if "END" not in line and "CONECT" not in line]

            for line in pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                    try:
                        resnum = line[22:27].strip()
                        atomnum = line[7:11].strip()
                        chain = line[21]
                        resname = line[17:20]
                        resnums.append(resnum)
                        atomnums.append(atomnum)
                        chains.append(chain)
                        resnames.append(resname)
                        if resname == "HOH":
                            water = f"{chain}:{resnum}"
                            if water not in self.user_waters and water not in self.water_to_exclude:
                                self.water_to_exclude.append(water)
                    #Line too short - Remarks pdb
                    except IndexError:
                        pass

        if self.n_waters < 1:
            return
        else:
            lig_length = resnames.count(self.ligand_residue)
            resnums = [int(num) for num in resnums if num]
            max_resnum = max(resnums)
            water_resnums = []
            water_chain = chains[0]  # water chain = 1st protein chain
            atomnum = max([int(num) for num in atomnums if num]) + 1 + lig_length
    
            water = cs.water * self.n_waters * n_inputs
            for inp in range(n_inputs):
                for n in range(self.n_waters):
                    O_coords = Vector([np.random.randint(0, 100) for i in range(3)])
                    H1_coords = O_coords + Vector(0.757, 0.586, 0.0)
                    H2_coords = O_coords + Vector(-0.757, 0.586, 0.0)
                    water_coords = water_coords + [list(O_coords)] + [list(H1_coords)] + [list(H2_coords)]
    
                    max_resnum += 1  # each water must have a different residue number
                    water_resnums = water_resnums + [max_resnum] * 3
                max_resnum += 1
    
            water_atomnums = [atomnum + j for j in range(self.n_waters * 3 * n_inputs)]
    
            # PDB lines - water
            water_output = []
            for atom, num, resnum, coord in zip(water, water_atomnums, water_resnums, water_coords):
                coord = ["{:7.4f}".format(c) for c in coord]
                coord = " ".join(coord)
                water_output.append(atom.format(num, water_chain, resnum, coord))
    
            sliced_water_output = []
            for i in range(0, len(water_output), self.n_waters * 3):
                sliced_water_output.append(water_output[i:i + self.n_waters * 3])
    
            # loop over minimisation inputs
            for inp, w in zip(self.input_pdbs, sliced_water_output):
                new_protein_file = inp
                # add water to PDB
                with open(inp, "w+") as file:
                    for line in pdb_lines:
                        file.write(line)
                    file.write("\n")
                    for line in w:
                        file.write(line)

                    file.write("END")
                # load again with Biopython
                parser = PDBParser()
                structure = parser.get_structure("complex", new_protein_file)
                water_list = []
                protein_list = Selection.unfold_entities(structure, "A")
                temp_protein_file = os.path.join(os.path.dirname(inp),
                                                 os.path.basename(inp).replace(".pdb", "_temp.pdb"))
                
                for res in structure.get_residues():
                    resnum = res._id[1]
                    if res.resname == 'HOH':
                        if resnum not in resnums:
                            water_list = water_list + Selection.unfold_entities(res, "A")
    
                # check for water contacts
                contacts5 = []
                for w in water_list:
                    contacts5 = contacts5 + NeighborSearch(protein_list).search(w.coord, 5.0, "A")
                contacts5 = [c for c in contacts5 if c not in water_list]  # exclude "self" contacts
    
                # translate water, if needed
                while contacts5:
                    contacts5 = []
                    for w_ in water_list:
                        x, y, z = w_.coord
                        w_.set_coord([x - 5, y, z])
                        contacts5 = contacts5 + NeighborSearch(protein_list).search(w_.coord, 5.0, "A")
                        contacts5 = [c for c in contacts5 if c not in water_list]
    
                # save final output
                io = PDBIO()
                io.set_structure(structure)
                io.save(temp_protein_file)
                output.append(new_protein_file)

                new_water_lines = []
                with open(temp_protein_file, "r") as temp:
                    temp_lines = temp.readlines()
                    for line in temp_lines:
                        if line[17:20].strip() == "HOH" and int(line[22:27].strip()) not in resnums:
                            line = line.replace(line[7:11], str(int(line[7:11]) + lig_length))
                            if line[12:15] == "2HW":
                                line = line.strip("\n") + "\nTER\n"
                            new_water_lines.append(line)

                del new_water_lines[-1] #Last biopython line is a not need it TER

                with open(new_protein_file, "w+") as file:
                    for line in pdb_lines:
                        file.write(line)
                    if not line.startswith("TER"):
                        file.write("TER\n")
                    for line in new_water_lines:
                        file.write(line)
                    for line in conect:
                        file.write(line)
                    file.write("\n")
                    file.write("END")
                os.remove(temp_protein_file)

def ligand_com(refinement_input, ligand_chain):
    parser = PDBParser()
    output = []
    refinement_input = glob.glob(refinement_input)

    for inp in refinement_input:
        structure = parser.get_structure("inp", inp)
        mass = 0.0
        com = np.zeros(3)
        for res in structure.get_residues():
            if res.resname == ligand_chain:
                for atom in res.get_atoms():
                    com = com + np.array(list(atom.get_vector())) * atom.mass
                    mass += atom.mass
                com = com / mass

        output.append(com.tolist())

    return output
