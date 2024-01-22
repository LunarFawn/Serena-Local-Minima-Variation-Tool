"""
File for managing analysis of dataset pnas.2112979119.sd01 from https://www.pnas.org/doi/full/10.1073/pnas.2112979119
"""

from enum import Enum
from typing import List
from pathlib import Path
from pandas import DataFrame

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.weighted_structures import WeightedStructure

from serena.interfaces.Sara2_API_Python3 import Sara2API, puzzleData
from serena.interfaces.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.interfaces.nupack4_0_28_wsl2_interface import (
    MaterialParameter, 
    NUPACK4Interface,
    NupackSettings
)

from serena.bin.backup_serena import ArchiveSecondaryStructureList



class ArchiveFlow(Enum):
    GET="GET"
    PUT="PUT"


class ProcessPNAS():
    
    def __init__(self) -> None:
        self._vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()
        self._nupack4:NUPACK4Interface = NUPACK4Interface()
    
    @property
    def vienna2_fmn_hack(self)->Vienna2FMNInterface:
        return self._vienna2_fmn_hack
    
    @property
    def  nupack4(self)->NUPACK4Interface:
        return self._nupack4
    
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

        
    def get_nupack_ensemble_structs(self, pnas_dataset_path:Path, round:str, sublab:str, nupack_settings:NupackSettings, do_archive:bool, archive_path:str=None)->List[Sara2StructureList]:
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        for design in puzzle_data.designsList:
            design_id= str(design.design_info.DesignID)
            sequence = design.design_info.Sequence
            fold_change = design.wetlab_results.FoldChange
            eterna_score = design.wetlab_results.Eterna_Score
            folding_subscore = design.wetlab_results.Folding_Subscore
            switch_subscore = design.wetlab_results.Switch_Subscore
            baseline_subscore = design.wetlab_results.Baseline_Subscore
            
            
            structs:Sara2StructureList = self.nupack4.get_subopt_energy_gap(material_param=nupack_settings.material_param,
                                                                    temp_C=nupack_settings.temp_C,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=nupack_settings.kcal_span_from_mfe,
                                                                    )
            
            self.archive_ensemble_structs(design_id=design_id,
                                          structs=[structs],
                                          dest_folder=archive_path,
                                          flow=ArchiveFlow.PUT)
    
    def archive_ensemble_structs(self,design_id:str, structs:List[Sara2StructureList], dest_folder:Path, flow:ArchiveFlow):
        pass
    
    def switchyness_analysis(self, designs_structures:List[Sara2StructureList]):
        pass
    

def get_nupack_ensemble_structs():
    pass

def archive_ensemble_structs():
    pass

def switchyness_analysis():
    pass