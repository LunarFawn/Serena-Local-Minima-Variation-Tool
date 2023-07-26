"""
File to hold the comparison structures files
"""

from dataclasses import dataclass
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.weighted_structures import WeightedStructureResult


@dataclass
class ComparisonResult():
    comp_struct: str = ''
    unbound_mfe_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()
    bound_mfe_struct: Sara2SecondaryStructure = Sara2SecondaryStructure()
    num_bound:float = -1
    num_unbound:float = -1
    num_both:float = -1
    num_dot:float = -1
    num_nucs:int = -1



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
        compared_struct:str = ''            

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
            
            compared_struct = compared_struct + comp_nuc_symbol
        
        compared_data: ComparisonResult = ComparisonResult(comp_struct=compared_struct,
                                                                           unbound_mfe_struct=unbound_struct,
                                                                           bound_mfe_struct=bound_struct,
                                                                           num_bound=num_bound,
                                                                           num_unbound=num_unbound,
                                                                           num_both=num_both,
                                                                           num_dot=num_dot,
                                                                           num_nucs=nuc_count)
        
        return compared_data

