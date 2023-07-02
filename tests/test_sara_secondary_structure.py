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
    stack_energgy:float - -20
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
    assert empty_secondary_structure.structure == -10

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
def secondary_structure_list():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    free_energy:float = -30
    stack_energgy:float - -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)


def test_default_new_secondary_struct_list(empty_secondary_structure_list:Sara2StructureList):
    assert empty_secondary_structure_list.mfe_structure == ''
    


































"""
Local Minima Structure Variation Data
Creation Date=2023-04-09 23:04:53.341105
---------------------------------------
***DESIGN INFO***
Design Name = Sara mod 1 of 6423414 #SaraFilterEverything-Good
DesignID = 6458932
Lab Name = Same State NG 1
Sequence = GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA
Eterna_Score = 100.0
FoldChange = 27.05
2nd State Target Structure = ........(((......(((.............))).....)))........................................
2nd State Folded Structure = ((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).
2nd State Folded Oligo Energy = -27.4
Energy Span from MFE = 7
Energy span units = .5
---------------------------------------
***RAW DATA***
Kcal,LMSV_U_mfe,LMSV_U_rel,LMSV_US_target,LMSV_US_folded
-31.679424285888672,1.0,1.0,48.0,22.0
-31.179424285888672,1.0,1.0,48.0,22.0
-30.679424285888672,2.0,0.0,51.0,23.0
-30.179424285888672,4.875,5.625,48.0,22.625
-29.679424285888672,4.571428571428571,4.857142857142858,47.0,22.71428571428571
-29.179424285888672,7.6875,8.0625,46.4375,17.625
-28.679424285888672,7.628571428571428,12.485714285714286,46.1142857142857,19.942857142857143
-28.179424285888672,8.92537313432836,19.95522388059703,46.67164179104479,20.014925373134336
-27.679424285888672,10.55882352941176,20.382352941176485,45.87254901960785,18.9607843137255
-27.179424285888672,10.880434782608695,19.266304347826082,45.88586956521739,19.30978260869564
-26.679424285888672,11.376811594202893,21.602898550724642,45.6608695652174,20.34202898550724
-26.179424285888672,12.073883161512022,22.18213058419245,45.687285223367695,20.096219931271488
-25.679424285888672,12.798024149286492,19.895718990120745,45.54774972557628,19.941822173435792
-25.179424285888672,13.47939262472885,18.55097613882864,45.09110629067244,19.798264642082437
---------------------------------------
***METRICS***
Polymorphicity Level (2kcal to end of sample) = 10.592338900810564
LMV_US_folded at folded with lignad/oligo energy = 19.30978260869564
---------------------------------------
EOF
"""