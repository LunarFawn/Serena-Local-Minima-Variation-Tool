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

from serena.structures import SingleEnsembleGroup, MultipleEnsembleGroups, Sara2SecondaryStructure, Sara2StructureList, EVResult, EV, LocalMinimaVariation
from serena.ensemble_variation import LMV_Token, LMV_ThreadProcessor, LMV_Shuttle


@dataclass
class WeightedStructureResult():
    ensemble_goup:SingleEnsembleGroup = SingleEnsembleGroup()
    weighted_struct: Sara2SecondaryStructure = Sara2SecondaryStructure()
    result_line: str = ''

@dataclass
class WeightedComparisonResult():
    comp_struct: str = ''
    unbound_mfe_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()
    bound_mfe_struct: Sara2SecondaryStructure = Sara2SecondaryStructure()
    num_bound:float = -1
    num_unbound:float = -1
    num_both:float = -1
    num_dot:float = -1

@dataclass
class WeightedValueRatios():
    unbound_to_total_nuc_ratio:float = -1
    last_unbound_to_current_unbound_ratio: float = -1
    bound_to_total_nuc_ratio:float = -1
    current_bound_to_last_bound_ratio: float = -1
    bound_to_unbound_ratio: float = -1


@dataclass
class WeightedLocalMinimaVariation():
    local_minima_variation_weighted_struct:EV = EV()
    local_minima_variation_unbound_struct:EV = EV()
    local_minima_variation_bound_struct:EV = EV()
    local_minima_variation_weighted_span_struct:EV = EV()

@dataclass
class WeightedScores():
    base_switch_score:float = 0
    penalties:int = 0
    bonuses:int = 0
    total_switch_score:int = 0
    is_on_switch:bool = False
    is_off_switch:bool = False
    is_functional_switch:bool = False
    is_low_foldchange:bool = False
    is_medium_foldchange:bool = False
    is_high_foldchange:bool = False

@dataclass
class WeightedGroupResult():
    weighted_struct_result: WeightedStructureResult = WeightedStructureResult()
    weighted_comparision: WeightedComparisonResult = WeightedComparisonResult()
    ratios: WeightedValueRatios = WeightedValueRatios()
    lmv: WeightedLocalMinimaVariation = WeightedLocalMinimaVariation()
    scores: WeightedScores = WeightedScores()

class EnsembleAndWeightedStructures():   

    def __init__(self) -> None:
        self._groups:List[SingleEnsembleGroup] = []
        self._comp_structs: List[WeightedStructureResult] = []
        self._num_groups:int = 0
    
    def add_group(self, group: SingleEnsembleGroup, comp_struct:WeightedStructureResult):
        self._groups.append(group)
        self._comp_structs.append(comp_struct)
        self._num_groups = self._num_groups + 1
    
    def get_group_and_comp_struct(self, index:int):
        return self._groups[index], self._comp_structs[index]
    
    @property
    def num_groups(self):
        return self._num_groups
    
    @property
    def comp_structs(self):
        return self._comp_structs

    @property
    def groups(self):
        return self._groups



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
                comp_struct = self.compair_weighted_structure(ensemble.multi_state_mfe[0], ensemble.multi_state_mfe[1], 
                                                                          weighted_struct, ensemble.group.nuc_count)                    
            else:
                comp_struct = "EMPTY GROUP"
        except Exception as error:
            comp_struct = f'BAD GROUP Error:{error}'
        
        result_line: str = f'{round(kcal_start,2)} to {round(kcal_stop,2)} kcal:   {comp_struct}'
        
        weighted_result: WeightedStructureResult = WeightedStructureResult(weighted_struct=weighted_struct, 
                                                         comp_struct=comp_struct, 
                                                         result_line=result_line)
        
        return weighted_result



    def make_weighted_struct(self, ensemble: SingleEnsembleGroup):
        structure_list: Sara2StructureList = ensemble.group
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

        saraw_weighted_struct: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=ensemble.group.sara_stuctures[0].sequence,
                                                                                    structure=weighted_structure,)
    
        weighted_result: WeightedStructureResult = WeightedStructureResult(weighted_struct=saraw_weighted_struct,
                                                                           ensemble_goup=ensemble)

        return weighted_result
    
    def compair_weighted_structure(self, unbound_mfe_struct:Sara2SecondaryStructure, bound_mfe_struct:Sara2SecondaryStructure, weighted_result:WeightedStructureResult, nuc_count:int):
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
            weighted_nuc:str = weighted_result.weighted_struct[nuc_index]
            unbound_nuc:str = unbound_mfe_struct.structure[nuc_index]
            bound_nuc: str = bound_mfe_struct.structure[nuc_index]

            comp_nuc_symbol:str = ''

            if weighted_nuc == bound_nuc and weighted_nuc != unbound_nuc:
                comp_nuc_symbol = bound
                num_bound += 1
            elif weighted_nuc != bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = unbound
                num_unbound += 1
            elif weighted_nuc == bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = both
                num_both += 1
            else:
                comp_nuc_symbol = dot
                num_dot += 1
            
            compared_struct = compared_struct + comp_nuc_symbol
        
        compared_data: WeightedComparisonResult = WeightedComparisonResult(comp_struct=compared_struct,
                                                                           unbound_mfe_struct=unbound_mfe_struct,
                                                                           bound_mfe_struct=bound_mfe_struct,
                                                                           num_bound=num_bound,
                                                                           num_unbound=num_unbound,
                                                                           num_both=num_both,
                                                                           num_dot=num_dot)
        
        return compared_data
    
    def process_ensemble_comp_lmv(weighted_groups_list:List[WeightedStructureResult], mfe_struct: Sara2SecondaryStructure, folded_struct: Sara2SecondaryStructure):

        lmv_results:List[WeightedLocalMinimaVariation] = []
        ensemble_groups: List[Sara2StructureList] =  []
        comp_structures: List[Sara2SecondaryStructure] = []
        rel_structs: List[Sara2SecondaryStructure] = []
  
        for weighted_group in weighted_groups_list:    
            ensemble_groups.append(weighted_group.ensemble_goup.group)
            comp_structures.append(weighted_group.weighted_struct)
            rel_structs.append(weighted_group.ensemble_goup.group.sara_stuctures[0].structure)

        #first run the MFE structs          

        lmv_mfe_thread: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=ensemble_groups,mfe_stuct=mfe_struct)
        lmv_mfe_thread_result:LMV_Token = lmv_mfe_thread.run_LMV()

        #now run comp structs
        lmv_comp_thread: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=ensemble_groups,comp_struct_list_option=comp_structures)
        lmv_comp_thread_result:LMV_Token = lmv_comp_thread.run_LMV()

        #now run rel structs
        lmv_rel_thread: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=ensemble_groups,mfe_struct=rel_structs)
        lmv_rel_thread_result:LMV_Token = lmv_rel_thread.run_LMV()

        #now run folded structs
        lmv_folded_thread: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=ensemble_groups,mfe_struct=rel_structs)
        lmv_folded_thread_result:LMV_Token = lmv_folded_thread.run_LMV()

        comp_ev_list_target: List[EV] = lmv_mfe_thread_result.group_results

        #print(comp_ev_list_target)

        ev_comp = comp_ev_list_target[0].ev_normalized
        ev_comp_limit: float = 25
        ev_mfe = group_ev_list_mfe[group_index].ev_normalized

        diff_limit:float = 1

        
        #if group_index == 1:
        if mfe_pronounced_first_group is True:
            diff:float = round(ev_mfe,2) - round(ev_comp,2)
            if round(ev_comp,2) < round(ev_mfe,2) and diff >= diff_limit:
                is_off_on_switch = True
                modifier = modifier + '+++'
                found_bound_list.append(group_index)
                if stop_diff is False:
                    found_bound_index = group_index
                stop_diff = True

        #if group_index == 0:
        diff = round(ev_comp,2) - round(ev_mfe,2)
        if round(ev_mfe,2) <= round(ev_comp,2):# and (diff >= diff_limit or diff == 0):
            mfe_pronounced_first_group = True

    def do_calculations_group(self, current_compared_data: WeightedComparisonResult, last_compared_data: WeightedComparisonResult, weighted_lmv:WeightedLocalMinimaVariation, raw_current_goup:SingleEnsembleGroup):
        start_group_mfe:float = raw_current_goup.kcal_start
        modifier= ''
        end_group_mfe:float = raw_current_goup.kcal_end
        folded_kcal:float = raw_current_goup.multi_state_mfe_kcal[1]
        bond_range_start:float = folded_kcal - 3
        bond_range_end:float = folded_kcal + 3
        last_unbound:float=last_compared_data.num_unbound
        last_bound:float=last_compared_data.num_bound
        is_functional_switch = False
        is_powerful_switch = False
        is_good_switch = False
        unbound_to_total_ratio:float = 0
        bound_ratio: float = 0
        last_unbound_ratio = 0
        last_bound_ratio = 0
        unbound = current_compared_data.num_unbound
        bound = current_compared_data.num_bound
        if unbound != 0:
            last_unbound_ratio = last_unbound/unbound 
            bound_ratio = bound/unbound
        if last_bound != 0:
            last_bound_ratio = bound/last_bound 
        unbound_to_total_ratio = unbound/raw_current_goup.group.nuc_count

        score:int = 0
        bonus:int = 0

        if start_group_mfe >= bond_range_start and start_group_mfe <= bond_range_end and end_group_mfe >= bond_range_start and end_group_mfe <= bond_range_end:
                    #if folded_kcal >=start_group_mfe and folded_kcal <= end_group_mfe:
                        is_in_bound_range = True
                        modifier = '***'

        bound_stats: str = f'BURatio:{round(bound_ratio,1)}, BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)} B:{bound}, U:{unbound}'
                    
        limit: float = 1.5 

        if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.25 and is_in_bound_range is True:
            is_good_switch = True
            score = score +1
        
        if last_unbound_ratio >= limit and last_bound_ratio >= limit and bound_ratio >=2 and is_in_bound_range is True:
            is_powerful_switch = True
            bonus = bonus +1

        if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and unbound_to_total_ratio <=.2 and is_in_bound_range is True:
            is_powerful_switch = True
            bonus = bonus +1

        if bound_ratio >=  limit and unbound_to_total_ratio <=.15 and is_in_bound_range is True:
            is_powerful_switch = True
            bonus = bonus +1

        total_score = score+ bonus
        weighted_scores: WeightedScores = WeightedScores(is_functional_switch=is_good_switch,
                                                         is_high_foldchange=is_powerful_switch,
                                                         base_switch_score=score,
                                                         bonuses=bonus,
                                                         total_switch_score=total_score
                                                         )
