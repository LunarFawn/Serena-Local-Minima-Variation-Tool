"""
Stuff I spun up just for eternacon presentation
"""


from dataclasses import dataclass
from serena.utilities.Sara2_API_Python3 import Sara2API, puzzleData

@dataclass
class PNASWetLab():
    pass

class run_eternacon():

    def __init__(self) -> None:
        pass

    def run(self):
        """
        run this
        """
        
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

        #should be able to for loop through all the designs now

