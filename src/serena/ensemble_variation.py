"""
This is the file to just get plain old ensemble variation for a list
of structures that make up an ensemble
"""


import attrs
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EnsembleVariation, EVResult


@attrs.define
class RunEnsembleVariation():

    def ev_from_list_strings_mfe(primary_structure:str, mfe_secondary_structure:str, secondary_structures_list: List[str])->float:
        #create mfe secondary structure
        mfe_structure:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=primary_structure,
                                                                        structure=mfe_secondary_structure)
        
        #populate secondary structure list
        structure_list:Sara2StructureList = Sara2StructureList()
        for structure in secondary_structures_list:
            structure_list.add_structure(Sara2SecondaryStructure(sequence=primary_structure,
                                                                 structure=structure))
            
        ensemble_variation:EnsembleVariation = EnsembleVariation()
        ev:EV = ensemble_variation.ensemble_variation_algorithm(kcal_group_structures_list=structure_list,
                                                        ref_structure=mfe_structure)
        
        return ev.ev_normalized