"""
Code for getting weighted structures.
"""

from typing import List, Dict
import struct
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
import threading
import time
import numpy as np
import collections

from serena.structures import SingleEnsembleGroup, MultipleEnsembleGroups, Sara2SecondaryStructure, Sara2StructureList, EVResult, WeightedStructureData

@dataclass
class WeightedResult():
    weighted_struct: str = ''
    comp_struct: str = ''
    result_line: str = ''

class WeightedStructures():
    """
    Class for struCts
    """

    def __init__(self) -> None:
        pass

    def process_ensemble_group(self, ensemble: SingleEnsembleGroup, kcal_start:float, kcal_stop:float, is_folded_mfe:bool, num_states: int):
        if num_states != 2:
            message :str =  "only supports 2 states right now"
            raise Exception(message)

        #for 2 states do this        
        group = ensemble.group
        comp_struct:str =''
        result_line:str = ''
        try:
            if group.num_structures > 0:
                weighted_struct = self.make_weighted_struct(group)
                comp_struct = self.compair_weighted_structure(ensemble.multi_state_mfe_struct[0], ensemble.multi_state_mfe_struct[1], 
                                                                          weighted_struct, ensemble.group.nuc_count)                    
            else:
                comp_struct = "EMPTY GROUP"
        except Exception as error:
            comp_struct = f'BAD GROUP Error:{error}'
        
        result_line: str = f'{round(kcal_start,2)} to {round(kcal_stop,2)} kcal:   {comp_struct}'
        
        weighted_result: WeightedResult = WeightedResult(weighted_struct=weighted_struct, 
                                                         comp_struct=comp_struct, 
                                                         result_line=result_line)
        
        return weighted_result



    def make_weighted_struct(self, structure_list: Sara2StructureList):
        is_bond_value: int = 2
        not_bond_value: int = -1

        nuc_poistion_values: List[int] = []
        nuc_pairs_comp_list: List[List[str]] = []
        good_nucs_each_pos: List[bool] = []

        struct_count: int = structure_list.num_structures

        for nucIndex in range(structure_list.nuc_count):
            nuc_poistion_values.append(0)
            pairs_list: List[str] = []            
            nuc_pairs_comp_list.append(pairs_list)
            #good_nucs_each_pos.append(False)

        for struct in structure_list.sara_stuctures:
            for nucIndex in range(structure_list.nuc_count):
                nuc_bond_type:str = struct.structure[nucIndex]
                nuc_pairs_comp_list[nucIndex].append(nuc_bond_type)
                adder: int = 0
                if nuc_bond_type == '.':
                    adder = not_bond_value
                else:
                    adder = is_bond_value
                nuc_poistion_values[nucIndex] = nuc_poistion_values[nucIndex] + adder
        
        #now record if the nuc position has a weghted bond
        for nucIndex in range(structure_list.nuc_count):
            is_weighted_bond=False
            if nuc_poistion_values[nucIndex] > struct_count:
                is_weighted_bond = True
            good_nucs_each_pos.append(is_weighted_bond)

        """
        for nucIndex in range(structure_list.nuc_count):
            nuc_value: float = float(nuc_poistion_values[nucIndex])
            #worked out this algotithm one night.. idk
            num_bonds_found:float = nuc_value / 2 + 1
            min_good_bonds: float = (num_bonds_found * 2) - ((nuc_value-num_bonds_found) * (-1))
            is_weighted_bond=False
            if num_bonds_found >= min_good_bonds and num_bonds_found > 0:
                is_weighted_bond = True
            good_nucs_each_pos.append(is_weighted_bond)
        """

        weighted_structure:str = ''
        for nucIndex in range(structure_list.nuc_count):
            is_bonded = good_nucs_each_pos[nucIndex]
            new_counter: collections.Counter = collections.Counter(nuc_pairs_comp_list[nucIndex])
            most_common_char: str= '.'
            if is_bonded is True:
                #most_common_char = '|'
                new_char:str = new_counter.most_common(2)[0][0]
                length = len(new_counter.most_common(2))
                if new_char == '.' and length > 1:
                    #then get second most common
                    new_char = new_counter.most_common(2)[1][0]
                most_common_char = new_char
            weighted_structure = weighted_structure + most_common_char
        
        return weighted_structure
    
    def compair_weighted_structure(self, unbound_mfe_struct:str, bound_mfe_struct:str, weighted_struct:str, nuc_count:int):
        """
        Compaire the weighted structure against the folded and not-folded mfe's.
        If a element is present in the folded mfe then it gets a '-'
        if element is in unbound only then it gets a '|'.
        The idea is that if you have a straight line in the list then it is very close to the
        folded mfe and if it is not straight then it is more like the unbound mfe.
        """
        unbound:str = '|'
        bound:str = '-'
        both:str = '+'
        dot:str = '.'
        compared_struct:str = ''            

        for nuc_index in range(nuc_count):
            weighted_nuc:str = weighted_struct[nuc_index]
            unbound_nuc:str = unbound_mfe_struct[nuc_index]
            bound_nuc: str = bound_mfe_struct[nuc_index]

            comp_nuc_symbol:str = dot

            if weighted_nuc == bound_nuc and weighted_nuc != unbound_nuc:
                comp_nuc_symbol = bound
            elif weighted_nuc != bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = unbound
            elif weighted_nuc == bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = both
            
            compared_struct = compared_struct + comp_nuc_symbol
        
        return compared_struct
    
    def calculate_performance_group(self, weighted_data: WeightedStructureData):
        group = weighted_data.raw_group
        folded_2nd_state_structure = weighted_data.bound_mfe_dot_paren_struct
        bound: int = 0
        unbound: int= 0
        print("weighted structs per group")
        mfe_energy = weighted_data.unbound_mfe_kcal
        start_group_mfe:float = mfe_energy + 0.5
        end_group_mfe:float = weighted_data.kcal_stop
        folded_kcal = weighted_data.bound_mfe_kcal
        bond_range_start:float = folded_kcal - 3
        bond_range_end:float = folded_kcal + 3
        last_unbound:float=0
        last_bound:float=0
        is_functional_switch = False
        is_powerful_switch = False
        is_good_switch = False
     
        comp_struct:str =''
        result:str = ''
        is_in_bound_range: bool = False
        modifier:str=''
        try:
            if group.num_structures > 0:
                new_struct = self.make_weighted_struct(group)
                comp_struct, bound, unbound = self.compair_weighted_structure(group.sara_stuctures[0].structure, folded_2nd_state_structure, new_struct, group.nuc_count)                    
                if start_group_mfe >= bond_range_start and start_group_mfe <= bond_range_end and end_group_mfe >= bond_range_start and end_group_mfe <= bond_range_end:
                #if folded_kcal >=start_group_mfe and folded_kcal <= end_group_mfe:
                    is_in_bound_range = True
                    modifier = '***'
            else:
                comp_struct = "no structures in kcal group"
        except Exception as error:
            comp_struct = f'bad list Error:{error}'
        unbound_to_total_ratio:float = 0
        bound_ratio: float = 0
        last_unbound_ratio = 0
        last_bound_ratio = 0
        if unbound != 0:
            last_unbound_ratio = last_unbound/unbound 
            bound_ratio = bound/unbound
        if last_bound != 0:
            last_bound_ratio = bound/last_bound 
        unbound_to_total_ratio = unbound/group.nuc_count

        bound_stats: str = f'BURatio:{round(bound_ratio,1)}, BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)} B:{bound}, U:{unbound}'
        last_unbound = unbound
        last_bound = bound
        line: str = f'{modifier} {round(start_group_mfe,2)} to {round(end_group_mfe,2)} kcal: {bound_stats}  {comp_struct}'
        print (line)

            
        limit: float = 1.5 

        if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.25 and is_in_bound_range is True:
            is_good_switch = True
        
        if last_unbound_ratio >= limit and last_bound_ratio >= limit and bound_ratio >=2 and is_in_bound_range is True:
            is_powerful_switch = True

        if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.2 and is_in_bound_range is True:
            is_powerful_switch = True

        if bound_ratio >=  limit and unbound_to_total_ratio <=.15 and is_in_bound_range is True:
            is_powerful_switch = True
            
        if is_good_switch is True or is_powerful_switch is True:
            print("Functional Switch")
            if is_powerful_switch is True:
                print('High Fold Change Predicted')
            else:
                print("Low fold change predicted")
        else:
            print("Bad Switch")
    
