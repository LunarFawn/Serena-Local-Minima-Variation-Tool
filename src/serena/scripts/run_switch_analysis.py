"""
Pyhton file that provides the classes neccessary to perform
local minima varitation calululations. This code is writtent to
be as agnostic to the source of data as possibel, but was developed using
nupack4, fyi.
copyright 2023 Jennifer Pearl
"""

from nupack import *
from dataclasses import dataclass
import math
import copy
from enum import Enum
#import matplotlib.pyplot as plt
#from matplotlib.ticker import StrMethodFormatter
#import matplotlib
from typing import List
from datetime import datetime
import numpy as np


from serena.apps.original_weighted_analysis import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult

debug:bool = False


from bisect import bisect_left

class SwitchPrediction(Enum):
    BAD = 0,
    FUNCTIONAL = 1
    EXCELLENT = 2,
    UNKOWN = 3,
    NONE = 4

@dataclass
class PredictionReponse():
    prediction:SwitchPrediction
    foldchange:float
    message:str
    raw_scores:List[float]
    num_structs:List[int]   

class OriginalSwitchAnalysis():

    def __init__(self) -> None:
        pass

    def take_closest(self, myList, myNumber):
        """
        Assumes myList is sorted. Returns closest value to myNumber.

        If two numbers are equally close, return the smallest number.
        """
        pos = bisect_left(myList, myNumber)
        if pos == 0:
            return myList[0]
        if pos == len(myList):
            return myList[-1]
        before = myList[pos - 1]
        after = myList[pos]
        if after - myNumber < myNumber - before:
            return after
        else:
            return before
        
    def do_switch_analysis(self, sequence, fmn_struct, fmn_struct_free_energy, span, units, manual:bool = False):
      
        target = '........(((......(((.............))).....)))........................................'
      
        if manual is True:
            print("Enter single strand RNA sequence")
            sequence = input()

            print("Enter target structure")
            target = '........(((......(((.............))).....)))........................................' #input()

            print("Enter predicted 2nd state folded structure")
            fmn_struct = input()

            print("Enter Energy of folded structure with ligand/oligo bound")
            fmn_struct_free_energy = float(input())

            print("Enter Kcal delta span to look at")        
            span = input()
            print(f'span is {span}')

            print("Enter kcal unit to plot by")
            units = '1' #input()
            print(f'units is {units}')

        EV_test: EnsembleVariation = EnsembleVariation()
        temp_list: List[int] = [36, 37, 38]
        score_list:List[float] = []
        raw_scores:List[float] = []
        num_structs:List[int] = []

        score: float = 0
        for temp in temp_list: 
            
            value, num_structs = EV_test.process_ensemble_variation(sequence, int(span), float(units), fmn_struct, target, fmn_struct_free_energy, temp)
            score = score + value
            score_list.append(value)
            raw_scores.append(value)
            num_structs.append(value)

        
        num_scores: int = len(temp_list)
        print(f'Raw score is {score} of {len(temp_list)}')
        
        num_zero_values: int = score_list.count(0)
        modified_score: float = (score - (num_zero_values)) / num_scores
        print(f'modified_score is {modified_score}')
        
        do_offset: bool = False
        if do_offset == True:
            max_fold_score:float = 3    
            if modified_score > (max_fold_score):
                #its probbaly only going to be in 2nd state state so give penatly
                offset = (modified_score - max_fold_score)
                modified_score = max_fold_score - offset

        predicted_foldchange_message:str = ''
        prediction: SwitchPrediction = SwitchPrediction.NONE
        why_not:float = (27/2.5)*modified_score
        predicted_foldchange:float = (27/3.5)*modified_score
        print(f'Predicted fold change is ~{predicted_foldchange}')
        #if modified_score > 3.5:
        #    predicted_foldchange = f'Good and Bad Switch SchrodenState Predicted. Fold performance Predicted to be so high (20s) that it might be too high. \nMight only form in 2nd state with low basescore and thus have bad fold change in wetlab. Probn KDOFF of ~150 and KDON off ~20ish or KDOFF of ~500 and KDON of ~25'
        #else:
        if predicted_foldchange > 15:
            predicted_foldchange_message = f'Good Switch Predicted. Fold change predicterd to be {predicted_foldchange} +/-6'
            prediction = SwitchPrediction.EXCELLENT
        elif predicted_foldchange == -1:
            predicted_foldchange_message = f'Unkown Switch Predicted. score of 0 for all temperatures so could be unkown fold to energy model or a really bad design'
            prediction = SwitchPrediction.UNKOWN
        elif predicted_foldchange < 1:
            predicted_foldchange_message = f'Bad Switch Predicted'
            prediction = SwitchPrediction.BAD
        else:
            predicted_foldchange_message = f'Underperforming Switch Predicted. Fold change predicterd to be {predicted_foldchange} +/-6'
            prediction = SwitchPrediction.FUNCTIONAL
        print(predicted_foldchange_message)

        response: PredictionReponse = PredictionReponse(prediction=prediction,
                                                        foldchange=predicted_foldchange,
                                                        message=predicted_foldchange_message,
                                                        raw_scores=raw_scores,
                                                        num_structs=num_structs)
        return response




