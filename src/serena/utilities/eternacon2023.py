"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
import statistics
import openpyxl
import pandas as pd
from pandas import DataFrame
from typing import List, Dict

from serena.utilities.Sara2_API_Python3 import Sara2API, puzzleData
from serena.utilities.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList, ComparisonStructures
from serena.utilities.weighted_structures import WeightedStructure
from serena.scripts.run_switch_analysis import OriginalSwitchAnalysis, PredictionReponse
from serena.utilities.logging import PNASAnalysisLogging

@dataclass
class PNASWetLab():
    pass

class run_eternacon():

    def __init__(self) -> None:
        pass

    def run_nupack_vienna(self):
        """
        run this
        """
        
        vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()

        pnas_path:str = ''
        pnas_round101_sheet:str = ''
        sublab_name:str = ''
        
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_path,
                                                      designRound_sheet=pnas_round101_sheet,
                                                      sublab_name=sublab_name)


        predicted_foldchange_list: Dict[int,float] = {}
        avg_raw_score_list: Dict[int,float] = {}
        avg_num_structures_list: Dict[int,float] = {}

        raw_score_36_list: Dict[int,float] = {}
        raw_score_37_list: Dict[int,float] = {}
        raw_score_38_list: Dict[int,float] = {}

        num_structures_36_list: Dict[int,float] = {}
        num_structures_37_list: Dict[int,float] = {}
        num_structures_38_list: Dict[int,float] = {}


        for design in puzzle_data.designsList:
            design_id= design.DesignInfo.DesignID
            sequence = design.DesignInfo.Sequence
            fold_change = design.wetlabResults.FoldChange
            eterna_score = design.wetlabResults.Eterna_Score
            folding_subscore = design.wetlabResults.Folding_Subscore
            switch_subscore = design.wetlabResults.Switch_Subscore
            baseline_subscore = design.wetlabResults.Baseline_Subscore

            #this is the fmn bound mfe struct, subopt list and weighted struck
            fmn_struct = vienna2_fmn_hack.rnafold_fmn(sequence)
            fmn_subopt = vienna2_fmn_hack.rnasubopt_fmn(sequence)
            fmn_weighted_struct: WeightedStructure = WeightedStructure(fmn_subopt)

            switch:OriginalSwitchAnalysis = OriginalSwitchAnalysis()

            analysis:PredictionReponse = switch.do_switch_analysis(sequence=sequence,
                                                        fmn_struct=fmn_struct,
                                                        fmn_struct_free_energy=0,
                                                        span=7,
                                                        units=1,
                                                        manual=False)
            
            predicted_foldchange_list[design_id] = analysis.foldchange
            
            
            avg_raw_score_list[design_id] = statistics.fmean(analysis.raw_scores)
            avg_num_structures_list[design_id] = statistics.fmean(analysis.num_structs)
            
            raw_score_36_list[design_id] = analysis.raw_scores[0]
            raw_score_37_list[design_id] = analysis.raw_scores[1]
            raw_score_38_list[design_id] = analysis.raw_scores[2]
            
            num_structures_36_list[design_id] = analysis.num_structs[0]
            num_structures_37_list[design_id] = analysis.num_structs[1]
            num_structures_38_list[design_id] = analysis.num_structs[2]
            
        pandas_sheet['PredictedFoldChange'] = predicted_foldchange_list
        
        pandas_sheet['AvgSwitchScore'] = avg_raw_score_list
        pandas_sheet['AvgNumSstruct'] = avg_num_structures_list
        
        pandas_sheet['36Deg_SwitchScore'] = raw_score_36_list
        pandas_sheet['37Deg_SwitchScore'] = raw_score_37_list
        pandas_sheet['38Deg_SwitchScore'] = raw_score_38_list

        pandas_sheet['36Deg_NumStructs'] = num_structures_36_list
        pandas_sheet['37Deg_NumStructs'] = num_structures_37_list
        pandas_sheet['38Deg_NumStructs'] = num_structures_38_list
     
        print("Its done!!!")

