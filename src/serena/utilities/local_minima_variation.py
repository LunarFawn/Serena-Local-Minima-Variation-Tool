"""
File to hold the local minima variation code
"""

from datetime import datetime
from typing import List, Dict
from dataclasses import dataclass

from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_variation import EV_Token, EV, EVResult
from serena.utilities.thread_manager import EV_ThreadProcessor
from serena.utilities.weighted_structures import WeightedEnsembleResult

@dataclass
class ComparisonLMV():
    lmv_comp:EV = EV()
    lmv_mfe:EV = EV()
    lmv_rel:EV = EV()

@dataclass
class ComparisonLMVResponse():
    lmv_comps:List[ComparisonLMV]
    

class LocalMinimaVariation():

    def __init__(self) -> None:
        pass

    def get_multi_group_lmv(self, ensemble: MultipleEnsembleGroups, reference_structure:Sara2SecondaryStructure):        
        print(f'Begining LMV processing at {datetime.now()}')
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,comp_structure=reference_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        lmv_results: EVResult = result_thread_LMV.ev_results
        return lmv_results

    def get_single_group_lmv(self, ensemble_group: SingleEnsembleGroup, reference_structure:Sara2SecondaryStructure):
        print(f'Begining LMV processing at {datetime.now()}')
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=[ensemble_group.group],
                                                              comp_structure=reference_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        lmv_results: EVResult = result_thread_LMV.ev_results
        return lmv_results
    
    

        
        

