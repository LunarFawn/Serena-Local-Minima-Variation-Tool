import pytest
from typing import List

from serena.generate_structures import SecondaryStructures
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList

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

def test_make_secondary_structure_list(secondary_structure_3:Sara2SecondaryStructure, secondary_structure_3_1: Sara2SecondaryStructure):
    new_struct_list:List[Sara2SecondaryStructure] = [secondary_structure_3, secondary_structure_3_1]
    generate_structs: SecondaryStructures = SecondaryStructures()
    secondary_structs_list: Sara2StructureList = generate_structs.make_secondary_strucuture_list(secondary_structures_list=new_struct_list)
    assert secondary_structs_list.num_structures == 2
    assert secondary_structs_list.sara_stuctures == new_struct_list

def test_make_single_ensemble_group():
    pass
