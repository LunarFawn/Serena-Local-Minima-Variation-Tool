"""
File to hold the comparison structures files
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure
#from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
#from serena.utilities.weighted_structures import WeightedEnsembleResult


@dataclass
class ComparisonNucCounts():
    bound_count:float = -1
    unbound_count:float = -1
    both_count:float = -1
    dot_count:float = -1
    num_nucs:int = -1

@dataclass
class ComparisonNucResults():
    comparison_nuc_counts: List[ComparisonNucCounts]

@dataclass
class ComparisonResult():
    #unbound_struct:Sara2SecondaryStructure 
    #bound_struct: Sara2SecondaryStructure
    #reference_struct: Sara2SecondaryStructure
    comp_struct: Sara2SecondaryStructure
    comp_counts: ComparisonNucCounts


class ComparisonStructures():

    def __init__(self) -> None:
        pass

    def compair_structures(self, unbound_struct:Sara2SecondaryStructure, bound_struct:Sara2SecondaryStructure, reference_struct:Sara2SecondaryStructure, nuc_count:int):
        """
        Compaire the weighted structure against the folded and not-folded mfe's.
        If a element is present in the folded mfe then it gets a '-'
        if element is in unbound only then it gets a '|'.
        The idea is that if you have a straight line in the list then it is very close to the
        folded mfe and if it is not straight then it is more like the unbound mfe.
        """
        unbound:str = '|'
        num_unbound:int = 0
        bound:str = '-'
        num_bound:int = 0
        both:str = '+'
        num_both:int = 0
        dot:str = '.'
        num_dot:int = 0
        temp_compared_struct:str = ''            

        for nuc_index in range(nuc_count):
            reference_nuc:str = reference_struct.structure[nuc_index]
            unbound_nuc:str = unbound_struct.structure[nuc_index]
            bound_nuc: str = bound_struct.structure[nuc_index]

            comp_nuc_symbol:str = ''

            if reference_nuc == bound_nuc and reference_nuc != unbound_nuc:
                comp_nuc_symbol = bound
                num_bound += 1
            elif reference_nuc != bound_nuc and reference_nuc == unbound_nuc:
                comp_nuc_symbol = unbound
                num_unbound += 1
            elif reference_nuc == bound_nuc and reference_nuc == unbound_nuc:
                comp_nuc_symbol = both
                num_both += 1
            else:
                comp_nuc_symbol = dot
                num_dot += 1
            
            temp_compared_struct = temp_compared_struct + comp_nuc_symbol
        
        comp_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=unbound_struct.sequence,
                                                                        structure=temp_compared_struct)
        
        comp_nuc_counts: ComparisonNucCounts = ComparisonNucCounts(bound_count=num_bound,
                                                                    unbound_count=num_unbound,
                                                                    both_count=num_both,
                                                                    dot_count=num_dot,
                                                                    num_nucs=nuc_count)
                
        compared_result: ComparisonResult = ComparisonResult(comp_struct=comp_struct,
                                                                #unbound_struct=unbound_struct,
                                                                #bound_struct=bound_struct,
                                                                #reference_struct=reference_struct,
                                                                comp_counts=comp_nuc_counts)
        
        return compared_result
    
    
    