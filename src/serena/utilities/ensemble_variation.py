"""
File for the ensemble variation code to live
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
from enum import Enum

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup

class SourceMFE(Enum):
    NONE = 0
    UNBOUND = 1
    BOUND = 2


@dataclass
class EV:
    ev_normalized: float = -1
    ev_ThresholdNorm: float = -1
    ev_structure: float = -1

@dataclass
class EVResult():
    #group_structs:SingleEnsembleGroup
    ev_values:List[EV]

@dataclass
class EVResult_old():
    groups_list : List[Sara2StructureList]
    groups_dict: Dict[int, Sara2StructureList]
    group_values: List[float]
    group_ev_list: List[EV]
    group_ev_dict: Dict[int,EV]


class EV_Token():
    def __init__(self, num_groups: int) -> None:
        self._group_results: List[EV] = num_groups * [EV()]
        self._group_dict: Dict[int,EV] = {}
        self._group_values: List[str] = num_groups * ['']
        self._group_done_status: List[bool] = num_groups * [False]
    
    @property
    def group_dict(self):
        return self._group_dict
        
    def set_group_dict(self, index:int, value:EV):
        self._group_dict[index]=value

    @property
    def group_results(self):
        return self._group_results

    @property
    def ev_results(self) -> EVResult:
        result: EVResult = EVResult(ev_values=self.group_results)
        return result
        
    def set_group_result(self, index:int, value:EV):
        self._group_results[index]=value
    
    @property
    def group_values(self):
        return self._group_values
        
    def set_group_values(self, index:int, value:str):
        self._group_values[index]=value

    @property
    def group_done_status(self):
        return self._group_done_status
        
    def set_group_done_status(self, index:int, state:bool):
        self._group_done_status[index]=state
    
    @property
    def is_done(self):
        is_completed:bool = False
        if self._group_done_status.count(False) == 0:
            #its done
            is_completed = True
        return is_completed


class EV_Shuttle():

    def __init__(self, structs_list: Sara2StructureList, mfe:Sara2SecondaryStructure, group_index:int, token:EV_Token) -> None:
        self._kcal_group_structures_list: Sara2StructureList = structs_list
        self._sara_mfestructure:Sara2SecondaryStructure = mfe
        self._group_index:int = group_index
        self._token:EV_Token = token
    
    @property
    def kcal_group_structures_list(self):
        return self._kcal_group_structures_list

    @kcal_group_structures_list.setter
    def kcal_group_structures_list(self, new_list: Sara2StructureList):
        self._kcal_group_structures_list = new_list
    
    @property
    def sara_mfestructure(self):
        return self._sara_mfestructure

    @sara_mfestructure.setter
    def sara_mfestructure(self, new_strucr: Sara2SecondaryStructure):
        self._sara_mfestructure = new_strucr
    
    @property
    def group_index(self):
        return self._group_index

    @group_index.setter
    def group_index(self, new_index: int):
        self._group_index = new_index
    
    @property
    def token(self):
        return self._token

    @token.setter
    def token(self, new_token: EV_Token):
        self._token = new_token

class EnsembleVariation():
    """
    Class that exposes the algorithm
    """

    def __init__(self) -> None:
        pass

    def thread_EV(self, shuttle: EV_Shuttle):
        total_EV_subscore1:int = 0
        structure_element_count = shuttle.kcal_group_structures_list.num_structures
        token:EV_Token = shuttle.token 
        group_num:int = shuttle.group_index
        if structure_element_count != 0:
            kcal_group_structures_list: Sara2StructureList = shuttle.kcal_group_structures_list

            sara_mfestructure:Sara2SecondaryStructure = shuttle.sara_mfestructure 
           
            
            #need to do each char abd then structure
            #walk through each nucleotide but first prep containers grab what is needed
            
            #setup constants
            nuc_count = kcal_group_structures_list.nuc_count
            structure_element_count = kcal_group_structures_list.num_structures

            ensembleVariation_score_temp = 0
            nucleotide_position_variation_basescores=[0]*nuc_count
            nucleotide_position_variation_subscores=[0]*nuc_count
            energydelta_individualVariationScore_list=[]
            
            #add the step to get nuc array here
            #get all the data out of it

            #first initialize the lists
            list_of_nuc_lists: List[List[str]] = []
    
            num_nucs: int = kcal_group_structures_list.nuc_count
            for index in range(num_nucs):
                temp_list:List[str] = []
                list_of_nuc_lists.append(temp_list)
                
            
            #now go throught everything
            for sara_structure in kcal_group_structures_list.sara_stuctures:
                for index in range(num_nucs):
                    character: str = sara_structure.structure[index]
                    list_of_nuc_lists[index].append(character)

            list_of_nuc_scores_base: List[int] = [0]*nuc_count
            list_of_nuc_scores_subscores: List[int] = [0]*nuc_count
            num_structs:int = kcal_group_structures_list.num_structures

            for nucIndex in range(nuc_count):
                mfe_nuc=sara_mfestructure.structure[nucIndex]
                num_chars = list_of_nuc_lists[nucIndex].count(mfe_nuc)
                num_diff:int = num_structs - num_chars
                list_of_nuc_scores_base[nucIndex] = num_diff
                list_of_nuc_scores_subscores[nucIndex] = list_of_nuc_scores_base[nucIndex] / structure_element_count 
            
            total_EV_subscore1 = sum(list_of_nuc_scores_subscores)
        else:
            total_EV_subscore1 = -1

        result: EV =  EV(ev_normalized=total_EV_subscore1, ev_ThresholdNorm=0, ev_structure=0)  
        token.group_results[group_num]= result
        token.group_dict[group_num] = result
        token.group_done_status[group_num] = True

