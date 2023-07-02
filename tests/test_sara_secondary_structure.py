import pytest
from typing import List

from serena.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        ComparisonStructures, 
                                        KcalRanges)


"""
Test sara2 secondary structure
"""

@pytest.fixture
def empty_secondary_structure():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2SecondaryStructure()

@pytest.fixture
def secondary_structure():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    free_energy:float = -30
    stack_energgy:float = -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)

def test_default_new_secondary_struct_list(empty_secondary_structure:Sara2SecondaryStructure):
    assert empty_secondary_structure.sequence == ''
    assert empty_secondary_structure.structure == ''
    assert empty_secondary_structure.freeEnergy == 0
    assert empty_secondary_structure.stackEnergy == 0
    assert empty_secondary_structure.nuc_count == 0

def test_setting_initial_secondary_struct(secondary_structure:Sara2SecondaryStructure):
    assert secondary_structure.sequence == 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    assert secondary_structure.structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert secondary_structure.freeEnergy == -30
    assert secondary_structure.stackEnergy == -20
    assert secondary_structure.nuc_count == 84

def test_setting_secondary_stuct_sequence(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.sequence = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    assert empty_secondary_structure.sequence == 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'

def test_setting_secondary_stuct_structure(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.structure = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert empty_secondary_structure.structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'

def test_setting_secondary_stuct_free_energy(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.freeEnergy = -10
    assert empty_secondary_structure.freeEnergy == -10

def test_setting_secondary_stuct_stack_energy(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.stackEnergy = -20
    assert empty_secondary_structure.stackEnergy == -20

def test_setting_secondary_stuct_nuc_count(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.sequence = 'GCC'
    assert empty_secondary_structure.nuc_count == 3



"""
Test sara2 secondary structure lists
"""

@pytest.fixture
def empty_secondary_structure_list():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2StructureList()

@pytest.fixture
def secondary_structure_1():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)

@pytest.fixture
def secondary_structure_2():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCG'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))))'
    free_energy:float = -50
    stack_energgy:float = -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)


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