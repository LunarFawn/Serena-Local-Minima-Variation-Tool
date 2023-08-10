from typing import List
import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonNucResults, ComparisonResult
from serena.utilities.weighted_structures import WeightedNucCounts,WeightedComparisonResult, WeightedStructure
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups

"""
Secondary Structure Fixtures
"""
@pytest.fixture
def kcal_range():
    return KcalRanges(start=1, stop=3)

@pytest.fixture
def kcal_range_2():
    return KcalRanges(start=2, stop=5)

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
def secondary_structures_list_2_item(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3)
    structure_list.add_structure(secondary_structure_4)
    return structure_list

@pytest.fixture
def secondary_structures_list_2_item_alt(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3_1)
    structure_list.add_structure(secondary_structure_5)
    return structure_list

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

"""
Fixtures for ensemble groups
"""

@pytest.fixture
def empty_single_ensemble_group():
    """
    Return a empty single ensemble group class
    """
    return SingleEnsembleGroup()

@pytest.fixture
def single_ensemble_group(secondary_structures_list_2_item:Sara2StructureList):
    """
    Return a empty single ensemble group class
    """
    ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
    ensemble_group.group = secondary_structures_list_2_item
    
    mfe_structs_list:List[str] = ['((..))','(...))']
    ensemble_group.multi_state_mfe_struct = mfe_structs_list
    
    mfe_kcal_list:List[float] = [-10,-20]
    ensemble_group.multi_state_mfe_kcal = mfe_kcal_list
    
    ensemble_group.kcal_end = 10
    ensemble_group.kcal_span = 20
    ensemble_group.kcal_start = 30
    return ensemble_group

@pytest.fixture
def single_ensemble_group_2(secondary_structures_list_2_item:Sara2StructureList):
    """
    Return a empty single ensemble group class
    """
    ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
    ensemble_group.group = secondary_structures_list_2_item
    
    mfe_structs_list:List[str] = ['(....)','..()..']
    ensemble_group.multi_state_mfe_struct = mfe_structs_list
    
    mfe_kcal_list:List[float] = [-30,-40]
    ensemble_group.multi_state_mfe_kcal = mfe_kcal_list
    
    ensemble_group.kcal_end = 40
    ensemble_group.kcal_span = 50
    ensemble_group.kcal_start = 60
    return ensemble_group



@pytest.fixture
def empty_multiple_ensemble_groups():
    """
    Return a empty multiple ensemble group class
    """
    return MultipleEnsembleGroups()

@pytest.fixture
def multiple_ensemble_groups(secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    """
    Returns a multiple ensemble groups class with
    values provided at instantiation
    """
    return MultipleEnsembleGroups(non_switch_kcal=10,
                                    non_switch_struct=secondary_structure_4,
                                    switched_kcal=20,
                                    switched_struct=secondary_structure_5)
