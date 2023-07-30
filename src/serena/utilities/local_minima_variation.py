"""
File to hold the local minima variation code
"""

from datetime import datetime
from typing import List, Dict
from dataclasses import dataclass

from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_variation import SourceMFE, EV_Token, EV, EVResult
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

@dataclass
class ReferenceStructures():
    mfe_structure:Sara2SecondaryStructure
    weighted_structures: WeightedEnsembleResult
    

class LocalMinimaVariation():

    def __init__(self) -> None:
        pass

    def get_multi_group_lmv(self, ensemble: MultipleEnsembleGroups, reference_structure:Sara2SecondaryStructure):        
        #now process all the groups
       
        #source_mfe: SourceMFE = SourceMFE.NONE

        #should be able to populate before hand and add to the group stuff i am working on
        #seams to push for the need a bit more
        print(f'Begining LMV processing at {datetime.now()}')
        
        #single_ensemble_group: List[Sara2StructureList] = [ensemble.group]
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                              comparison_structure=reference_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        group_ev_list: List[EV] = result_thread_LMV.group_results
        group_ev_dict: Dict[int,EV] = result_thread_LMV.group_dict

        #need to change to make this list part of the call
        #lmv_results: List[EVResult] = []
        lmv_results: EVResult = result_thread_LMV.ev_results
        #for group_index in range(len(ensemble.groups)):
        #    group: SingleEnsembleGroup = ensemble.groups[group_index]
        #    ev:EV = group_ev_list[group_index]
        #    ev_result:EVResult = EVResult(group_structs=group,
        #                                  ev_values=ev)
        #    lmv_results.append(ev_result)

        return lmv_results


    def get_single_group_lmv(self, ensemble_group: SingleEnsembleGroup, reference_structure:Sara2SecondaryStructure):
        print(f'Begining LMV processing at {datetime.now()}')
        #single_ensemble_group: List[Sara2StructureList] = [ensemble.group]
        LMV_Thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=[ensemble_group.group],
                                                              comparison_structure=reference_structure)
        result_thread_LMV:EV_Token = LMV_Thread.run_EV()
        lmv_results: EVResult = result_thread_LMV.ev_results
        return lmv_results
    
    def get_relative_mutli_group_lmv(self, ensemble: MultipleEnsembleGroups):
        #result_list:List[EVResult]= []
        ev_values:List[EV] = []
        for group in ensemble.groups:
            ref_structure:Sara2SecondaryStructure = group.group.mfe_structure
            ev_result:EVResult = self.get_single_group_lmv(ensemble_group=group,
                                                        reference_structure=ref_structure)
            ev_values.append(ev_result.ev_values[0])
        result: EVResult = EVResult(ev_values=ev_values)
        return result
    
    def get_weighted_multi_group_lmv(self, ensemble: MultipleEnsembleGroups, weighted_structures: WeightedEnsembleResult):
        #result_list:List[EVResult]= []
        ev_values:List[EV] = []
        for group_index in range(len(ensemble.groups)):

            ref_structure:Sara2SecondaryStructure = weighted_structures.structs[group_index]
            ev_result:EVResult = self.get_single_group_lmv(ensemble_group=ensemble.groups[group_index],
                                                        reference_structure=ref_structure)
            ev_values.append(ev_result.ev_values[0])
        result: EVResult = EVResult(ev_values=ev_values)
        return result
    
    def process_serena_lmvs(self, ensemble: MultipleEnsembleGroups, ref_structures:ReferenceStructures):
        
        
        #first get mfe lmv then weighted for groups
        mfe_result:EVResult = self.get_multi_group_lmv(ensemble=ensemble,
                                                        reference_structure=ref_structures.mfe_structure)
        #now get ref ev
        rel_result:EVResult = self.get_relative_mutli_group_lmv(ensemble=ensemble)

        #now get weightedEV
        weight_result:EVResult = self.get_weighted_multi_group_lmv(ensemble=ensemble, 
                                                                   weighted_structures=ref_structures.weighted_structures)
        comparisons_lmv_response: List[ComparisonLMV] = []
        for group_index in range(len(ensemble.groups)):
            lmv_data:ComparisonLMV = ComparisonLMV()
            lmv_data.lmv_comp = weight_result[group_index]
            lmv_data.lmv_mfe = mfe_result[group_index]
            lmv_data.lmv_rel = rel_result[group_index]
            comparisons_lmv_response.append(lmv_data)
        
        serena_lmv_respone: ComparisonLMVResponse = ComparisonLMVResponse(lmv_comps=comparisons_lmv_response)
        return serena_lmv_respone

        
        

