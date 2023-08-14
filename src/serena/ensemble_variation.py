"""
This is the file to just get plain old ensemble variation for a list
of structures that make up an ensemble
"""


import attrs
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EnsembleVariation
from serena.utilities.weighted_structures import WeightedEnsembleResult



class RunEnsembleVariation():
    
    def get_ev_from_structures_list(self, structures_list:Sara2StructureList, mfe_structure:Sara2SecondaryStructure):
        ensemble_variation:EnsembleVariation = EnsembleVariation()
        ev:EV = ensemble_variation.ensemble_variation_algorithm(kcal_group_structures_list=structures_list,
                                                        ref_structure=mfe_structure)
        
        return ev.ev_normalized