# Serena-Local-Minima-Variation-Tool

Installation Instructions:
Serna is a pip installable package

To install package:
pip install .

To uninstall package pass:
pip uninstall serena-local-minima-variation-tool

User Instructions:

Serena consists of my take and thoughts on how to design the framework for transporting RNA bioinformatics to and from other applications, as well as provide access to novel RNA analysis algorithms and methodologies.

To call the framework:

import Serena

The Sara2 Secondary Structure is the backbone of the framework. It is an object that contains everything needed to represent a RNA secondary structure in dot bracket notation with the primary structure sequence and both free engery and stack energy.

To create a sara 2 secondary structure example from unit test:

from serena.generate_structures import SecondaryStructures

def test_make_secondary_structure():
    generate_structs: SecondaryStructures = SecondaryStructures()
    sequence:str = "ACGUAC"
    structure:str = '((()))'
    free_energy:float = -31
    stack_energy:float = -41
    secondary_structure:Sara2SecondaryStructure = generate_structs.make_secondary_structure(primary_structure=sequence,
                                                                    secondary_structure=structure,
                                                                    free_energy=free_energy,
                                                                    stack_free_energy=stack_energy)
    assert secondary_structure.structure == structure
    assert secondary_structure.sequence == sequence
    assert secondary_structure.nuc_count == 6
    assert secondary_structure.freeEnergy == free_energy
    assert secondary_structure.stackEnergy == stack_energy

The sara2 secondary structure list is the primary vehicle for performing analysis of the structural ensemble. It is populated with each secondary structure found during for the ensemble.

To create a sara 2 secondary structure list example from unit test:

from serena.generate_structures import SecondaryStructures

def test_make_secondary_structure_list(secondary_structure_3:Sara2SecondaryStructure, secondary_structure_3_1: Sara2SecondaryStructure):
    new_struct_list:List[Sara2SecondaryStructure] = [secondary_structure_3, secondary_structure_3_1]
    generate_structs: SecondaryStructures = SecondaryStructures()
    secondary_structs_list: Sara2StructureList = generate_structs.make_secondary_strucuture_list(secondary_structures_list=new_struct_list)
    assert secondary_structs_list.num_structures == 2
    assert secondary_structs_list.sara_stuctures == new_struct_list