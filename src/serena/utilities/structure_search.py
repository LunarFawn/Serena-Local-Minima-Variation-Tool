"""
File for handeling the searching for various types of structures
"""

from enum import Enum
from dataclasses import dataclass
from typing import List

import re
from re import Match, Pattern

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

@dataclass
class MoleculareSnareDef():
    first_half_molecule:str
    length_first_half:int
    first_molecule_match:Match
    second_half_molecule:str
    lenght_second_half:int
    second_molecule_match:Match
    static_stem_start_index:int
    static_stem_end_index:int
    static_stem_structure:str
    

@dataclass
class SnareResults():
    snare_dectected:bool
    snare_list:List[MoleculareSnareDef]
    number_snares:int

class MolecularSnare():
    
    def __init__(self) -> None:
        self.nupack:NUPACK4Interface = NUPACK4Interface()
    
    def find_two_segments_with_stem_snare(self, moleculte_binding_sequence:str, unbound_secondary_structure:Sara2SecondaryStructure, bound_secondary_structure:Sara2SecondaryStructure)->SnareResults:
        unbound_pairs_list:List[str] = self.nupack.sara2_pairs_list(secondary_structure=unbound_secondary_structure)
        bound_pairs_list:List[str] = self.nupack.sara2_pairs_list(secondary_structure=bound_secondary_structure)
        
        first_half_molecule:str = ''
        lenght_first_half:int = -1
        second_half_molecule:str = ''
        lenght_second_half:int = -1
        
        snare_dectected:bool = False
        snare_list:List[MoleculareSnareDef] = []
        number_snares:int = 0
        
        search_sequence:str = unbound_secondary_structure.sequence
        
            
        for index in reversed(range(len(moleculte_binding_sequence))):
            
            #only need to look at unbound as the primary structure will not change and
            #we only want to look at the unbound state representation for now
            
            #will start with frist half being bigger than second half
            #work your way up the sequence and check as you go
            first_half_molecule = moleculte_binding_sequence[:index]
            lenght_first_half = len(first_half_molecule)
            second_half_molecule = moleculte_binding_sequence[index:]
            lenght_second_half = len(second_half_molecule)
            
            first_half_moelcule_pattern:Pattern = re.compile(first_half_molecule)
            first_half_molecule_matches:List[Match] = first_half_moelcule_pattern.search(search_sequence)
            
            if first_half_molecule_bound_structure == None:
                #there are no matches so move on to the next one
                continue
            
            for first_half_match in first_half_molecule_matches:
                
                #only allow the search to be after the end of the first half
                second_half_molecule_pattern:Pattern = re.compile(first_half_molecule)
                second_half_molecule_matches:List[Match] = second_half_molecule_pattern.search(search_sequence, first_half_match.end())
                
                first_half_end_index:int = first_half_match.end()
                if second_half_molecule_matches == None:
                    # there are no matches so skip to the next match
                    continue
                
                for second_half_match in second_half_molecule_matches:
                    
                    # if second_half_match.start() < first_half_end_index:
                    #     #the second half is before the first half so should not form
                    #     #could form if psuedoknot, but very unlikly and not the focus of
                    #     #this search for now
                    #     continue
                    
                    second_half_start_index:int = second_half_match.start()
                    
                    bound_structure_segment:str = bound_secondary_structure.structure[first_half_end_index:second_half_start_index]
                    unbound_structure_segment:str = unbound_secondary_structure.structure[first_half_end_index:second_half_start_index]
                    #do +1 and -1 due to the ends being able to be paired as part of the stacks that are associated with the snare
                    #this makes it so that we are only checking if the snare is part of a loop for verification of existence
                    first_half_molecule_bound_structure:str = bound_secondary_structure.structure[first_half_match.start()+1:first_half_match.end()-1]
                    second_half_molecule_bound_structure:str = bound_secondary_structure.structure[second_half_match.start()+1:second_half_match.end()-1]
                    
                    if unbound_structure_segment == bound_structure_segment:
                        if first_half_molecule_bound_structure.count('.') == len(first_half_molecule_bound_structure) and second_half_molecule_bound_structure.count('.') == len(second_half_molecule_bound_structure):
                            #it is most likely a snare, but now need to do a double check and make sure the stem is a static stem
                            #and is forming the stack for the snare and not bound somewhere else. maybe still keep track of all the posible configurations you can think there
                            #could be
                            
                            
                            
                            snare_dectected = True
                            number_snares += 1  
                            # snare_list.append(unbound_structure_segment)
                            new_snare:MoleculareSnareDef = MoleculareSnareDef(first_half_molecule=first_half_molecule,
                                                                            first_molecule_match=first_half_match,
                                                                            second_half_molecule=second_half_molecule,
                                                                            second_molecule_match=second_half_match,
                                                                            static_stem_start_index=first_half_end_index,
                                                                            static_stem_end_index=second_half_start_index,
                                                                            static_stem_structure=unbound_structure_segment,
                                                                            length_first_half=lenght_first_half,
                                                                            lenght_second_half=lenght_second_half) 
                            snare_list.append(new_snare)        
                # else:
                #     #there is no second half found in the snare
                #     pass
                        
        #now return result
        snare_search_results:SnareResults = SnareResults(snare_dectected=snare_dectected,
                                                         snare_list=snare_list,
                                                         number_snares=number_snares)
        
        return snare_search_results    
            
            
            # second_half_matches:List[Match] = re.search(second_half_molecule, search_sequence)
            
            # for second_half_match_index in reversed(range(len(second_half_matches))):
                
            
            #     first_half_sequence:str = search_sequence[:second_half_matches[second_half_match_index].start()]
            
            #     # second_half_sequence:str = search_sequence[start_index_second_half_sequence:]
            #     first_half_matches:List[Match] = re.search(first_half_molecule, first_half_sequence)
            
            #     for match in first_half_matches:
                    
            
            
            
            
            
            
                         
            # if second_half_molecule in search_sequence: # first_half_molecule in search_sequence and
                
            #     #start with the second half as this is bigger and prone to les
            #     start_index_second_half_sequence:int =  search_sequence.rfind(second_half_molecule)
                
            #     number_second_half_molecule:int = search_sequence.count(second_half_molecule)
                
            #     first_half_sequence:str = search_sequence[:start_index_second_half_sequence]
                
            #     # second_half_sequence:str = search_sequence[start_index_second_half_sequence:]
            #     # first_half_matches:List[Match] = re.search(first_half_molecule, first_half_sequence)
                
            #     number_first_half_molecule_in_first_half_sequence:int = first_half_sequence.count(first_half_molecule)
                
            #     if number_first_half_molecule_in_first_half_sequence < 1:
            #         #dont do this then so contirnue
            #         continue
                
            #     for first_index in range(number_first_half_molecule_in_first_half_sequence):
            #         #get the start index of the last occurance of the first half molecule
            #         start_first_half_molecule:int = first_half_sequence.rfind(first_half_molecule)
                    
            #         end_index_first_half_molecule:int = start_first_half_molecule + lenght_first_half
            #         start_index_second_half_molecule:int = start_index_second_half_sequence

            #         bound_investigation_structure:str = bound_secondary_structure.structure[end_index_first_half_molecule:start_index_second_half_molecule]
            #         unbound_investigation_structure:str = unbound_secondary_structure.structure[end_index_first_half_molecule:start_index_second_half_molecule]
                
            #         if bound_secondary_structure == unbound_secondary_structure:
            #             #we have found a molecular snare
            #             snare_dectected = True
                        
            #             #how many are there now
            #             number_snares += 1
            #         else:
            #             #if it is not found then go to the next occurance of the first half
            #             #make the new first half sequence stop where the last one started.
            #             #we are working our way down
            #             first_half_sequence = first_half_sequence[:start_first_half_molecule]
                        

                
            
        
        
        