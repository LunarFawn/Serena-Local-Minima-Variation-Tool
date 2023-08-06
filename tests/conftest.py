import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonNucResults, ComparisonResult

@pytest.fixture()
def empty_secondary_structure():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2SecondaryStructure()


@pytest.fixture()
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

@pytest.fixture()
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

@pytest.fixture
def empty_secondary_structure_list():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2StructureList()

@pytest.fixture
def empty_kcal_range():
    """
    Returns an initialized kcal range with default values
    """
    return KcalRanges()

@pytest.fixture
def empty_comparison_nuc_count():
    """
    Returns an empty comparions nuc pair dataclass
    """
    return ComparisonNucCounts()

@pytest.fixture
def comparison_nuc_count():
    """
    Returns populated comparions nuc pair dataclass
    """
    return ComparisonNucCounts(bound_count = 1,
                                unbound_count=2,
                                both_count=3,
                                dot_count=4,
                                num_nucs=5)