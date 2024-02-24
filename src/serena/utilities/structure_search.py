"""
File for handeling the searching for various types of structures
"""

from enum import Enum
from dataclasses import dataclass
from typing import List, Dict

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
    snare_stem_nuc_count:int
    is_prime_stem:bool
    is_hairpin_stem:bool
    

@dataclass
class SnareResults():
    snare_dectected:bool
    snare_list:List[MoleculareSnareDef]
    number_snares:int

class MolecularSnare():
    
    def __init__(self) -> None:
        self.nupack:NUPACK4Interface = NUPACK4Interface()
    
    # def find_snare_binding_rotation(self, moleculte_binding_sequence:str, unbound_secondary_structure:Sara2SecondaryStructure, bound_secondary_structure:Sara2SecondaryStructure):
    #     pass    
            
    # def find_hairpin_moleculare_snare(self):
    #     pass
    
    def find_prime_moleculare_snare(self, moleculte_binding_sequence:str, unbound_secondary_structure:Sara2SecondaryStructure, bound_secondary_structure:Sara2SecondaryStructure)->SnareResults:
        """
        prime moleculare snare has the 5'3' static stem as the static stem for the snare. 
        
        There is also rotation of the molecule binding pattern. This results in the joint for the binding site
        being at one of the two ends of the snare loop 
        
        both sides of the bindin site need to bnd to forma loop that holds teh molecule. the 2nd sate will be a full loop but the 1st state
        by nature must not be a full loop and only one side will be static. that is the snare side. 
        """
        unbound_pairs_list:List[int] = self.nupack.sara2_pairs_list(secondary_structure=unbound_secondary_structure)
        bound_pairs_list:List[int] = self.nupack.sara2_pairs_list(secondary_structure=bound_secondary_structure)
        
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
                
                first_half_start_index:int = first_half_match.start()
                first_half_end_index:int = first_half_match.end()
                
                #only allow the search to be after the end of the first half
                second_half_molecule_pattern:Pattern = re.compile(second_half_molecule)
                second_half_molecule_matches:List[Match] = second_half_molecule_pattern.search(search_sequence, first_half_end_index)
                
                
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
                    second_half_end_index:int = second_half_match.end()
                    
                    bound_structure_segment_prime_side_1:str = bound_secondary_structure.structure[:first_half_start_index]
                    bound_structure_segment_prime_side_2:str = bound_secondary_structure.structure[second_half_end_index:]
                    unbound_structure_segment_prime_side_1:str = unbound_secondary_structure.structure[:first_half_start_index]
                    unbound_structure_segment_prime_side_2:str = unbound_secondary_structure.structure[second_half_end_index:]
                    
                    bound_structure_segment:str = bound_secondary_structure.structure[first_half_end_index:second_half_start_index]
                    unbound_structure_segment:str = unbound_secondary_structure.structure[first_half_end_index:second_half_start_index]
                    #do +1 and -1 due to the ends being able to be paired as part of the stacks that are associated with the snare
                    #this makes it so that we are only checking if the snare is part of a loop for verification of existence
                    first_half_molecule_bound_structure:str = bound_secondary_structure.structure[first_half_start_index+1:first_half_end_index-1]
                    second_half_molecule_bound_structure:str = bound_secondary_structure.structure[second_half_start_index+1:second_half_end_index-1]
                    
                    #first check that the molecule binding sites are only a loop
                    molecule_is_loop:bool = False
                    static_hairpin: bool = False
                    static_prime:bool = False
                    if first_half_molecule_bound_structure.count('.') == len(first_half_molecule_bound_structure) and second_half_molecule_bound_structure.count('.') == len(second_half_molecule_bound_structure):
                        #it is most likely a snare, but now need to do a double check and make sure the stem is a static stem
                        #and is forming the stack for the snare and not bound somewhere else. maybe still keep track of all the posible configurations you can think there
                        #could be
                        
                        #if it is a loop then now test for static stems on both sides
                        molecule_is_loop = True
                    
                    #now we need to test
                    side_A_nuc_pair_list:List[tuple] = []
                    side_B_nuc_pair_list:List[tuple] = []
                    
                    snare_stem_nuc_count:int = 0
                    if molecule_is_loop is True and unbound_structure_segment == bound_structure_segment:    
                        #this means that it has a non-prime snare or what I am calling a hairpin snare
                        static_hairpin = True
                        in_stem:bool = False
                        in_loop:bool = True
                        for side_A_nuc_index in range(second_half_start_index, second_half_start_index-1, -1):
                            side_A_bound_nuc_index_pair:int = bound_pairs_list[side_A_nuc_index]
                            side_A_unbound_nuc_index_pair:int = unbound_pairs_list[side_A_nuc_index]
                            
                            if side_A_unbound_nuc_index_pair == side_A_bound_nuc_index_pair:
                                #it is for sure static and the pairing is the same
                                new_tuple:tuple = (side_A_nuc_index, side_A_bound_nuc_index_pair)
                                side_A_nuc_pair_list.append(new_tuple)
                        
                        for side_B_nuc_index in range(first_half_end_index, first_half_end_index+1, 1):
                            side_B_bound_nuc_index_pair:int = bound_pairs_list[side_B_nuc_index]
                            side_B_unbound_nuc_index_pair:int = unbound_pairs_list[side_B_nuc_index]
                            
                            if side_B_unbound_nuc_index_pair == side_B_bound_nuc_index_pair:
                                #it is for sure static and the pairing is the same
                                new_tuple:tuple = (side_B_nuc_index, side_B_bound_nuc_index_pair)
                                side_B_nuc_pair_list.append(new_tuple)
                        
                        snare_stem_nuc_count = len(bound_structure_segment)
                                    
                        
                            # #check if it is in a loop and if so then continue walking the struct until hit a pair
                            # #this should be a reasonable distance from the molecule binding site
                            # if in_stem is True:
                            #     #it has been detected that we are now in the stem
                            #     pass
                            
                            # if bound_nuc_index_pair == nuc_index and unbound_nuc_index_pair == nuc_index:
                            #     in_loop = True
                            #     # this means that we are in a loop
                            #     if in_stem is True:
                            #         #this means that we have hit a break in the stem
                            #     else:
                                    
                                    
                    if molecule_is_loop is True and bound_structure_segment_prime_side_1 == unbound_structure_segment_prime_side_1 and bound_structure_segment_prime_side_2 == unbound_structure_segment_prime_side_2:
                        static_prime = True
                        for side_A_nuc_index in range(first_half_start_index, first_half_start_index-1, -1):
                            side_A_bound_nuc_index_pair:int = bound_pairs_list[side_A_nuc_index]
                            side_A_unbound_nuc_index_pair:int = unbound_pairs_list[side_A_nuc_index]
                            
                            if side_A_unbound_nuc_index_pair == side_A_bound_nuc_index_pair:
                                #it is for sure static and the pairing is the same
                                new_tuple:tuple = (side_A_nuc_index, side_A_bound_nuc_index_pair)
                                side_A_nuc_pair_list.append(new_tuple)
                        
                        for side_B_nuc_index in range(second_half_end_index, second_half_end_index+1, 1):
                            side_B_bound_nuc_index_pair:int = bound_pairs_list[side_B_nuc_index]
                            side_B_unbound_nuc_index_pair:int = unbound_pairs_list[side_B_nuc_index]
                            
                            if side_B_unbound_nuc_index_pair == side_B_bound_nuc_index_pair:
                                #it is for sure static and the pairing is the same
                                new_tuple:tuple = (side_B_nuc_index, side_B_bound_nuc_index_pair)
                                side_B_nuc_pair_list.append(new_tuple)
                        
                        snare_stem_nuc_count = len(bound_structure_segment_prime_side_1) + len(bound_structure_segment_prime_side_2)
                                            
                    if molecule_is_loop is True and (static_hairpin is True or static_prime is True):     
                        
                         #now compare the dictionaries. 
                        for side_nuc_index in range(len(side_A_nuc_pair_list)):
                            side_A_nuc:int = side_A_nuc_pair_list[side_nuc_index][0]
                            side_A_pair:int = side_A_nuc_pair_list[side_nuc_index][1]
                            side_B_nuc:int = side_B_nuc_pair_list[side_nuc_index][0]
                            side_B_pair:int = side_B_nuc_pair_list[side_nuc_index][1]
                            
                            
                                #this should be the either the start of the bindig ring
                                #or it is the end
                            if side_A_nuc == side_A_pair and side_B_nuc == side_B_pair:
                                #this first one is part of the loop so we can continue
                                if side_nuc_index == 0:
                                    continue
                                else:
                                    break
                            elif side_A_nuc == side_B_pair and side_A_pair == side_B_nuc:
                                snare_dectected = True
                                break
                            # elif side_nuc_index > 0:
                                
                            #     if 
                                    
                        if snare_dectected == True:                        
                            #if this all it true then we now need to 
                            #investigate the pairs list 3.14159265358979323846
                            
                            number_snares += 1  
                            # snare_list.append(unbound_structure_segment)
                            new_snare:MoleculareSnareDef = MoleculareSnareDef(first_half_molecule=first_half_molecule,
                                                                            first_molecule_match=first_half_match,
                                                                            second_half_molecule=second_half_molecule,
                                                                            second_molecule_match=second_half_match,                                                                            
                                                                            length_first_half=lenght_first_half,
                                                                            lenght_second_half=lenght_second_half,
                                                                            snare_stem_nuc_count=snare_stem_nuc_count,
                                                                            is_prime_stem=static_prime,
                                                                            is_hairpin_stem=static_hairpin) 
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
                        

                
            
        
        
        