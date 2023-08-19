"""
This file is the entry point to get lmv results for an ensemble
"""

from datetime import datetime
from typing import List, Dict
from dataclasses import dataclass


from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EnsembleVariation, EVResult, EV_Shuttle, EV_Token
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
from serena.utilities.thread_manager import EV_ThreadProcessor
from serena.utilities.weighted_structures import WeightedEnsembleResult
from serena.utilities.local_minima_variation import LocalMinimaVariation


class RunLocalMinimaVariation(LocalMinimaVariation):
   
    def get_relative_mutli_group_lmv(self, ensemble: MultipleEnsembleGroups):
        ev_values:List[EV] = []
        for group in ensemble.groups:
            ref_structure:Sara2SecondaryStructure = group.group.sara_stuctures[0]
            ev_result:EVResult = self.get_single_group_lmv(ensemble_group=group,
                                                        reference_structure=ref_structure)
            ev_values.append(ev_result.ev_values[0])
        result: EVResult = EVResult(ev_values=ev_values)
        return result
    
    def get_relative_multi_group_lmv_nupack():
        pass
    
    def get_weighted_multi_group_lmv(self, ensemble: MultipleEnsembleGroups, weighted_structures: WeightedEnsembleResult):
        ev_values:List[EV] = []
        for group_index in range(len(ensemble.groups)):

            ref_structure:Sara2SecondaryStructure = weighted_structures.structs[group_index]
            ev_result:EVResult = self.get_single_group_lmv(ensemble_group=ensemble.groups[group_index],
                                                        reference_structure=ref_structure)
            ev_values.append(ev_result.ev_values[0])
        result: EVResult = EVResult(ev_values=ev_values)
        return result
    
    def get_mfe_mult_group_lmv(self, ensemble: MultipleEnsembleGroups):
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                            comp_structure=ensemble.non_switch_state_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        lmv_results: EVResult = result_thread_LMV.ev_results
        return lmv_results
    
    