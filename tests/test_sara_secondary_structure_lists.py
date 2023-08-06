
import pytest

from serena.utilities. ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        )


"""
Test sara2 secondary structure lists
"""


def test_default_new_secondary_struct_list(empty_secondary_structure_list:Sara2StructureList):
    assert empty_secondary_structure_list.mfe_structure == ''
    assert empty_secondary_structure_list.mfe_freeEnergy == 0
    assert empty_secondary_structure_list.mfe_stackEnergy == 0
    assert empty_secondary_structure_list.nuc_count == 0
    assert empty_secondary_structure_list.sara_stuctures == []
    assert empty_secondary_structure_list.max_free_energy == 0
    assert empty_secondary_structure_list.min_free_energy == 0
    assert empty_secondary_structure_list.max_stack_energy == 0
    assert empty_secondary_structure_list.min_stack_energy == 0
    assert empty_secondary_structure_list.num_structures == 0 
    assert empty_secondary_structure_list.freeEnergy_span == 0
    assert empty_secondary_structure_list.stackEnergy_span == 0
    assert empty_secondary_structure_list.weighted_structure == ''

def test_add_sara_struct_sara_list(empty_secondary_structure_list:Sara2StructureList, secondary_structure_1: Sara2SecondaryStructure, secondary_structure_2: Sara2SecondaryStructure):
    empty_secondary_structure_list.add_structure(secondary_structure_1)
    empty_secondary_structure_list.add_structure(secondary_structure_2)
    assert empty_secondary_structure_list.sara_stuctures[0] == secondary_structure_1
    assert empty_secondary_structure_list.sara_stuctures[1] == secondary_structure_2
    assert empty_secondary_structure_list.mfe_structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert empty_secondary_structure_list.mfe_freeEnergy == -30
    assert empty_secondary_structure_list.mfe_stackEnergy == -10
    assert empty_secondary_structure_list.nuc_count == 84
    assert empty_secondary_structure_list.max_free_energy == -30
    assert empty_secondary_structure_list.min_free_energy == -50
    assert empty_secondary_structure_list.max_stack_energy == -10
    assert empty_secondary_structure_list.min_stack_energy == -20
    assert empty_secondary_structure_list.num_structures == 2 
    assert empty_secondary_structure_list.freeEnergy_span == 20
    assert empty_secondary_structure_list.stackEnergy_span == 10

def test_set_weighted_struct_sara_list(empty_secondary_structure_list:Sara2StructureList):
    empty_secondary_structure_list.weighted_structure = '.....'
    assert empty_secondary_structure_list.weighted_structure == '.....'