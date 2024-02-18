"""
File for handeling the searching for various types of structures
"""

from enum import Enum
from dataclasses import dataclass
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure

from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface

@dataclass
class PrimeNucCounts():
    prime_static_nuc_count_total: int = 0
    five_prime_nuc_count: int = 0
    three_prime_nuc_count: int = 0
    static_stem_nuc_count: int = 0
    five_prime_dot_count:int = 0
    three_prime_dot_count:int = 0
    static_loop_nuc_count:int = 0
    

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
        
        five_prime_dot_count:int = 0
        three_prime_dot_count:int = 0
        
        
        in_stack:bool = False
        for five_prime_index in range(num_nucs):
            unbound_char: str = unbound_structure.structure[five_prime_index]
            bound_char: str = bound_structure.structure[five_prime_index]
            if in_stack == True:
                if unbound_char == bound_char:
                    # if in_stack is True:
                    if unbound_char == ')' or unbound_char == '.':
                        break
                else:
                    break
                        
                    
            if unbound_char == '(' and bound_char == '(':
                in_stack = True
                five_prime_nuc_count += 1
            elif unbound_char == '.' or bound_char == '.':
                five_prime_dot_count +=1                    
                # else:
                #     break
        
        in_stack = False
        for three_prime_index in reversed(range(num_nucs)):
            unbound_char: str = unbound_structure.structure[three_prime_index]
            bound_char: str = bound_structure.structure[three_prime_index]
            if in_stack is True:    
                if unbound_char == bound_char:                
                    if unbound_char == '(' or unbound_char == '.':
                        break
                else:
                    break
                    
                
            if unbound_char == ')' and bound_char == ')':
                in_stack = True
                three_prime_nuc_count += 1
            elif unbound_char == '.' or bound_char == '.':
                three_prime_dot_count +=1
                        
                
            
        
        prime_static_nuc_count_total = five_prime_nuc_count + three_prime_nuc_count + three_prime_dot_count + five_prime_dot_count
        
        prime_count.prime_static_nuc_count_total = prime_static_nuc_count_total
        prime_count.five_prime_nuc_count = five_prime_nuc_count
        prime_count.three_prime_nuc_count = three_prime_nuc_count
        prime_count.static_stem_nuc_count = five_prime_nuc_count + three_prime_nuc_count #min([five_prime_nuc_count,three_prime_nuc_count]) * 2
        prime_count.three_prime_dot_count = three_prime_dot_count
        prime_count.five_prime_dot_count = five_prime_dot_count
        prime_count.static_loop_nuc_count = three_prime_dot_count + five_prime_dot_count
        # prime_count.static_loop_plus_nuc_count = max([five_prime_nuc_count,three_prime_nuc_count]) - min([five_prime_nuc_count,three_prime_nuc_count])
        return prime_count


class MolecularSnare():
    
    def __init__(self) -> None:
        self.nupack:NUPACK4Interface = NUPACK4Interface()
    
    def find_two_segments_with_stem_snare(self, moleculte_binding_sequence:str, unbound_secondary_structure:Sara2SecondaryStructure, bound_secondary_structure:Sara2SecondaryStructure):
        unbound_pairs_list:List[str] = self.nupack.sara2_pairs_list(secondary_structure=unbound_secondary_structure)
        bound_pairs_list:List[str] = self.nupack.sara2_pairs_list(secondary_structure=bound_secondary_structure)
        
        first_half_molecule:str = ''
        lenght_first_half:int = -1
        second_half_molecule:str = ''
        lenght_second_half:int = -1
        
        snare_dectected:bool = False
        
        search_sequence:str = unbound_secondary_structure.sequence
        
        for index in range(len(moleculte_binding_sequence)):
            
            #only need to look at unbound as the primary structure will not change and
            #we only want to look at the unbound state representation for now
            first_half_molecule = moleculte_binding_sequence[:index]
            lenght_first_half = len(first_half_molecule)
            second_half_molecule = moleculte_binding_sequence[index:]
            lenght_second_half = len(second_half_molecule)
                                    
            if second_half_molecule in search_sequence: # first_half_molecule in search_sequence and
                
                #start with the second half as this is bigger and prone to les
                start_index_second_half_sequence:int =  search_sequence.rfind(second_half_molecule)
                
                number_second_half_molecule:int = search_sequence.count(second_half_molecule)
                
                first_half_sequence:str = search_sequence[:start_index_second_half_sequence]
                
                # second_half_sequence:str = search_sequence[start_index_second_half_sequence:]
                number_first_half_molecule_in_first_half_sequence:int = first_half_sequence.count(first_half_molecule)
                
                if number_first_half_molecule_in_first_half_sequence < 1:
                    #dont do this then so contirnue
                    continue
                
                for first_index in range(number_first_half_molecule_in_first_half_sequence):
                    #get the start index of the last occurance of the first half molecule
                    start_first_half_molecule:int = first_half_sequence.rfind(first_half_molecule)
                    
                
                
                

                
            
        
        
        