"""
This is the file to just get plain old ensemble variation for a list
of structures that make up an ensemble
"""


from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EnsembleVariation
from serena.utilities.weighted_structures import WeightedEnsembleResult
from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, NupackSettings, MaterialParameter


class RunEnsembleVariation(EnsembleVariation):
    """
    Class to get the ensemble variation
    """
    def ev_from_structures_list(self, structures_list:Sara2StructureList, mfe_structure:Sara2SecondaryStructure):
        ev:EV = self.ensemble_variation_algorithm(kcal_group_structures_list=structures_list,
                                                        ref_structure=mfe_structure)
        return ev.ev_normalized
    
    def ev_from_nupack4(self, material:MaterialParameter, temp_C:int, span_from_mfe:int, sequence:str):
        nupack4: NUPACK4Interface = NUPACK4Interface()   
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material,
                                  temp_C=temp_C,
                                  sequence_string=sequence,
                                  energy_delta_from_MFE=span_from_mfe,
                                  )
        ensemble_variation:float = self.ev_from_structures_list(structures_list=structs, mfe_structure=structs.sara_stuctures[0])
        return ensemble_variation
        