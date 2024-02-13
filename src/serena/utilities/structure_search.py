"""
File for handeling the searching for various types of structures
"""

from enum import Enum
from dataclasses import dataclass

from serena.utilities.ensemble_structures import Sara2SecondaryStructure

@dataclass
class PrimeNucCounts():
    prime_static_nuc_count_total: int = 0
    five_prime_nuc_count: int = 0
    three_prime_nuc_count: int = 0
    static_stem_nuc_count: int = 0
    static_loop_plus_nuc_count:int = 0
    

class StaticSystemDetector():

    def __init__(self) -> None:
        pass
    
    def find_3prime_5prime_static_system(self, unbound_structure: Sara2SecondaryStructure, bound_structure:Sara2SecondaryStructure)->PrimeNucCounts:
        dot_bracket:str = unbound_structure.structure
        num_nucs:int = unbound_structure.nuc_count
        
        prime_count: PrimeNucCounts = PrimeNucCounts()
        prime_static_nuc_count_total: int = 0
        five_prime_nuc_count: int = 0
        three_prime_nuc_count: int = 0
        
        for five_prime_index in range(num_nucs):
            if unbound_structure.structure[five_prime_index] == bound_structure.structure[five_prime_index]:
                five_prime_nuc_count += 1
            else:
                break
        
        for three_prime_index in reversed(range(num_nucs)):
            if unbound_structure.structure[three_prime_index] == bound_structure.structure[three_prime_index]:
                three_prime_nuc_count += 1
            else:
                break
        
        prime_static_nuc_count_total = five_prime_nuc_count + three_prime_nuc_count
        
        prime_count.prime_static_nuc_count_total = prime_static_nuc_count_total
        prime_count.five_prime_nuc_count = five_prime_nuc_count
        prime_count.three_prime_nuc_count = three_prime_nuc_count
        prime_count.static_stem_nuc_count = min([five_prime_nuc_count,three_prime_nuc_count]) * 2
        prime_count.static_loop_plus_nuc_count = max([five_prime_nuc_count,three_prime_nuc_count]) - min([five_prime_nuc_count,three_prime_nuc_count])
        return prime_count

