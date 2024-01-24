"""
File for managing analysis of dataset pnas.2112979119.sd01 from https://www.pnas.org/doi/full/10.1073/pnas.2112979119
"""

from enum import Enum
from typing import List, Union
from pathlib import Path
from pandas import DataFrame
import os
from dataclasses import dataclass

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.weighted_structures import WeightedStructure

from serena.interfaces.Sara2_API_Python3 import Sara2API, puzzleData, DesignPerformanceData
from serena.interfaces.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.interfaces.nupack4_0_28_wsl2_interface import (
    MaterialParameter, 
    NUPACK4Interface,
    NupackSettings
)

from serena.bin.backup_serena_v2 import ArchiveSecondaryStructureList

from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, MaterialParameter, NupackSettings, EnsembleSwitchStateMFEStructs
from serena.analysis.ensemble_analysis import InvestigateEnsemble, InvestigateEnsembleResults
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from serena.utilities.logging_serena import PNASAnalysisLogging


class ArchiveFlow(Enum):
    GET="GET"
    PUT="PUT"

@dataclass
class ArchiveData():
    design_info:DesignPerformanceData = None
    nupack_settings:NupackSettings = None
    structs:Sara2StructureList = None
    fmn_folded_mfe:Sara2SecondaryStructure = None
    fmn_folded_weighted:Sara2SecondaryStructure = None
    
class ProcessPNAS():
    
    def __init__(self) -> None:
        self._vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()
        self._nupack4:NUPACK4Interface = NUPACK4Interface()
        self._scoreing:InvestigateEnsemble = InvestigateEnsemble()
    
    @property
    def vienna2_fmn_hack(self)->Vienna2FMNInterface:
        return self._vienna2_fmn_hack
    
    @property
    def nupack4(self)->NUPACK4Interface:
        return self._nupack4

    @property
    def scoring(self)->InvestigateEnsemble:
        return self._scoreing
    
    def get_fmn_state_fold(self, sequence:str, do_weighted:bool)->Sara2SecondaryStructure:
        
        #  = False
        struct_to_use:Sara2SecondaryStructure = ''
        if do_weighted is True:
        #this is the fmn bound mfe struct, subopt list and weighted struck
            fmn_subopt = self.vienna2_fmn_hack.rnasubopt_fmn(input_sequence=sequence,
                                                        do_fmn=True)
            fmn_weighted_struct: WeightedStructure = WeightedStructure()
            struct_to_use= fmn_weighted_struct.make_weighted_struct(fmn_subopt)
        else:
            struct_to_use = self.vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
                                                        do_fmn=True)
        return struct_to_use

        
    def record_nupack_ensemble_structs(self, pnas_dataset_path:Path, round:str, sublab:str, nupack_settings:NupackSettings, do_archive:bool, archive_path:str=None)->List[Sara2StructureList]:
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        pnas_data:List[Sara2StructureList] = []
        for design in puzzle_data.designsList:
            # design_id= str(design.design_info.DesignID)
            sequence = design.design_info.Sequence
            # fold_change = design.wetlab_results.FoldChange
            # eterna_score = design.wetlab_results.Eterna_Score
            # folding_subscore = design.wetlab_results.Folding_Subscore
            # switch_subscore = design.wetlab_results.Switch_Subscore
            # baseline_subscore = design.wetlab_results.Baseline_Subscore
            
            
            structs:Sara2StructureList = self.nupack4.get_subopt_energy_gap(material_param=nupack_settings.material_param,
                                                                    temp_C=nupack_settings.temp_C,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=nupack_settings.kcal_span_from_mfe,
                                                                    )
            pnas_data.append(structs)
            fmn_mfe:Sara2SecondaryStructure = self.get_fmn_state_fold(sequence=sequence,
                                                                      do_weighted=False)
            fmn_weighted_struct:Sara2SecondaryStructure = self.get_fmn_state_fold(sequence=sequence,
                                                                                  do_weighted=True)
            
            # if do_archive is True:
            if os.path.isdir(archive_path) is False:
                raise FileExistsError(f'File {archive_path} is not a valid path')
            archive_data:ArchiveData = ArchiveData(design_info=design,
                                                    nupack_settings=nupack_settings,
                                                    structs=structs,
                                                    fmn_folded_mfe=fmn_mfe,
                                                    fmn_folded_weighted=fmn_weighted_struct)
            
            self.archive_ensemble_data(dest_folder=archive_path,
                                        flow=ArchiveFlow.PUT,data=archive_data)

    def retrieve_archive_and_run_analysis(self, save_folder:Path, run_name:str, pnas_dataset_path:Path, round:str, sublab:str, do_weighted:bool, is_agressive:bool = False, max_num_structs:int= 500000, archive_path:str=None):
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        # pnas_data:List[Sara2StructureList] = []
        for design in puzzle_data.designsList:
            if os.path.isdir(archive_path) is False:
                raise FileExistsError(f'File {archive_path} is not a valid path')
            
            found_data:ArchiveData = self.archive_ensemble_data(dest_folder=archive_path,
                                        flow=ArchiveFlow.GET)
            
            if found_data.structs.num_structures < max_num_structs:
                struct_to_use:Sara2SecondaryStructure = found_data.fmn_folded_mfe
                if do_weighted is True:
                    struct_to_use = found_data.fmn_folded_weighted
                    
                reference_structures:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(switched_mfe_struct=struct_to_use,
                                                                                                non_switch_mfe_struct=found_data.structs.sara_stuctures[0])
                
                ensemble_groups: MultipleEnsembleGroups = self.nupack4.load_nupack_subopt_as_ensemble(span_structures=found_data.structs,
                                                                                                kcal_span_from_mfe=found_data.nupack_settings.kcal_span_from_mfe,
                                                                                                Kcal_unit_increments=found_data.nupack_settings.Kcal_unit_increments,
                                                                                                switch_state=reference_structures
                                                                                                )
                
                investigation_results:InvestigateEnsembleResults = self.scoring.investigate_and_score_ensemble(ensemble=ensemble_groups,
                                                                                                           is_aggressive=is_agressive)
                
                total_scores: float = 0
                #if investigation_results.basic_scores.total_score > 0:
                total_scores: float = investigation_results.basic_scores.total_score + investigation_results.advanced_scores.total_score
                #else:
                #    total_scores = 0 - investigation_results.advanced_scores.excess_struct_penalty
                
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore'] = total_scores
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore_NoExcessStructs'] = total_scores + investigation_results.advanced_scores.excess_struct_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Basic_score'] = investigation_results.basic_scores.total_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Functional_score'] = investigation_results.basic_scores.functional_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Powerfull_score'] = investigation_results.basic_scores.powerful_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Bonuses'] = investigation_results.basic_scores.bonuses
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'On_off_score'] = investigation_results.basic_scores.on_off_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Penalties'] = investigation_results.basic_scores.penalties
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = investigation_results.advanced_scores.total_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_bonus'] = investigation_results.advanced_scores.comp_bonus
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_penalty'] = investigation_results.advanced_scores.comp_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_bonus'] = investigation_results.advanced_scores.lmv_bonus
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_penalty'] = investigation_results.advanced_scores.lmv_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = investigation_results.advanced_scores.excess_struct_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'NumStructs'] =  investigation_results.number_structures
                
                
            else:
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore_NoExcessStructs'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Basic_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Functional_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Powerfull_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Bonuses'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'On_off_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Penalties'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_bonus'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_bonus'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = 100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'NumStructs'] =  found_data.structs.num_structures
            
            design_data_df:DataFrame = pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID]
            logging: PNASAnalysisLogging = PNASAnalysisLogging()
            if os.path.isdir(save_folder.as_posix()) == False:
                os.makedirs(save_folder.as_posix())
            logging.save_excel_sheet(design_data_df, save_folder, run_name)
            
    
    def archive_ensemble_data(self, dest_folder:Path, flow:ArchiveFlow, data: ArchiveData=None)->Union[None, ArchiveData]:
        
        backup_records:ArchiveSecondaryStructureList = ArchiveSecondaryStructureList(working_folder=dest_folder,
                                             var_name=str(data.design_info.design_info.DesignID),
                                             use_db=True)
        
           
        if flow == ArchiveFlow.PUT:
            backup_records.archive.archive_data = data
            return None
        elif flow == ArchiveFlow.GET:
            retrieved_archive:ArchiveData = backup_records.archive.archive_data 
            return retrieved_archive
                                 

    
    # def switchyness_analysis(self, designs_structures:List[Sara2StructureList]):
    #     pass
    

def get_nupack_ensemble_structs():
    pass

def archive_ensemble_structs():
    pass

def switchyness_analysis():
    pass