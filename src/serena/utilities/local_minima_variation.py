"""
File to hold the local minima variation code
"""

from datetime import datetime
from typing import List, Dict

from serena.utilities.ensemble_groups import MultipleEnsembleGroups
from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_variation import SourceMFE, EV_Token, EV, EVResult
from serena.utilities.thread_manager import EV_ThreadProcessor

class LocalMinimaVariation():

    def __init__(self) -> None:
        pass

    def get_local_minima_variation(self, ensemble: MultipleEnsembleGroups, comparison_structure:Sara2SecondaryStructure):        
        #now process all the groups
       
        #source_mfe: SourceMFE = SourceMFE.NONE

        #should be able to populate before hand and add to the group stuff i am working on
        #seams to push for the need a bit more
        print(f'Begining LMV processing at {datetime.now()}')
        
        #single_ensemble_group: List[Sara2StructureList] = [ensemble.group]
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                              comparison_structure=comparison_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        group_ev_list: List[EV] = result_thread_LMV.group_results
        group_ev_dict: Dict[int,EV] = result_thread_LMV.group_dict

        result_LMV: EVResult = EVResult(groups_list=ensemble.raw_groups, groups_dict=ensemble.groups_dict, 
                                              group_values=ensemble.group_values, group_ev_list=group_ev_list, 
                                              group_ev_dict=group_ev_dict)
        return result_LMV