"""
Main entry for the algorithm to determine the swithchyness 
of a sequence baed on analysis of the ensemble
"""

from serena.analysis.ensemble_analysis import InvestigateEnsemble, InvestigateEnsembleResults
from serena.interfaces.nupack4_0_28_wsl2_interface import MaterialParameter, NUPACK4Interface
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs

class RunInvestigateEnsemble(InvestigateEnsemble):
    """
    Class that is the main entry point for ensemble investigation
    """
    
    def investigate_and_score_ensemble_nupack(self,sequence:str, folded_referenec_struct:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, kcal_unit_increments: float = 1, aggressive:bool= False)->InvestigateEnsembleResults:
        """
        Use the nupack folding enginer to generate a MultipleEnsembleGroups from a sequence 
        and refence folded structure (folded mfe maybe?) and analyze it for switchyness score
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=structs.sara_stuctures[0],
                                                                                    switched_mfe_struct=Sara2SecondaryStructure(sequence=sequence,
                                                                                                                                structure=folded_referenec_struct))
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        results:InvestigateEnsembleResults = self.investigate_and_score_ensemble(ensemble=ensemble,
                                                                                 is_aggressive=aggressive)
        return results