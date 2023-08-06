from typing import List
import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonNucResults, ComparisonResult
from serena.utilities.weighted_structures import WeightedNucCounts,WeightedComparisonResult, WeightedStructure

"""
Secondary Structure Fixtures
"""
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

@pytest.fixture()
def secondary_structure_3():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '((.)))'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)

@pytest.fixture()
def secondary_structure_3_1():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '((..))'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)

@pytest.fixture()
def secondary_structure_4():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '..().)'
    free_energy:float = -50
    stack_energgy:float = -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   freeEnergy=free_energy,
                                   stackEnergy=stack_energgy)

@pytest.fixture()
def secondary_structure_5():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '(...))'
    free_energy:float = -40
    stack_energgy:float = -30
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
def secondary_structures_list(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3)
    structure_list.add_structure(secondary_structure_3_1)
    structure_list.add_structure(secondary_structure_4)
    structure_list.add_structure(secondary_structure_5)
    return structure_list

@pytest.fixture
def empty_kcal_range():
    """
    Returns an initialized kcal range with default values
    """
    return KcalRanges()



"""
Comparison Structure Fixtures
"""

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

@pytest.fixture
def comparison_nuc_count_2():
    """
    Returns populated comparions nuc pair dataclass
    """
    return ComparisonNucCounts(bound_count = 2,
                                unbound_count=3,
                                both_count=4,
                                dot_count=5,
                                num_nucs=5)                            

@pytest.fixture
def comparison_nuc_result(comparison_nuc_count:ComparisonNucCounts, comparison_nuc_count_2:ComparisonNucCounts):
    """
    Returns a comparison nuc result
    """
    comp_list:List[ComparisonNucResults] = []
    comp_list.append(comparison_nuc_count)
    comp_list.append(comparison_nuc_count_2)
    return ComparisonNucResults(comparison_nuc_counts=comp_list)

@pytest.fixture
def comparison_result(secondary_structure_3:Sara2SecondaryStructure, comparison_nuc_count:ComparisonNucCounts):
    """
    Returns a comparison result (not the list one)
    """
    return ComparisonResult(comp_struct=secondary_structure_3,
                            comp_counts=comparison_nuc_count)


"""
Weighted Structure Fixtures
"""

@pytest.fixture
def empty_weighted_nuc_count():
    """
    Returns a empty weighted nuc count
    """
    return WeightedNucCounts()

@pytest.fixture
def weighted_nuc_count():
    """
    Returns a weighted nuc count populated at initialization
    """
    return WeightedNucCounts(num_unbound=1,
                            num_both=2,
                            num_bound=3,
                            num_dot=4,
                            num_nucs=5)

@pytest.fixture
def empty_weighted_comparison_result():
    """
    Return an empty comparison result
    """
    return WeightedComparisonResult()

@pytest.fixture
def weighted_struct_class():
    return WeightedStructure()