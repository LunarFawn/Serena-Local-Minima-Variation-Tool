
import pytest
from typing import List

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure, ComparisonStructures
from serena.utilities.weighted_structures import WeightedStructure


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

def test_new_weighted_structure():
    pass