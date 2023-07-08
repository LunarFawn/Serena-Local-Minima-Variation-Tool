"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
import openpyxl

from serena.utilities.Sara2_API_Python3 import Sara2API, puzzleData
from serena.utilities.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList, ComparisonStructures,
from serena.utilities.weighted_structures import WeightedStructure
from serena.scripts.run_switch_analysis import OriginalSwitchAnalysis, PredictionReponse

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
        
        new_sara:Sara2API = Sara2API()

        puzzle_data: puzzleData = new_sara.ProcessLab()

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

            wb = openpyxl.Workbook()
            ws_write = wb.create_sheet(0)
            ws_write.append(desing_header)
        #should be able to for loop through all the designs now

