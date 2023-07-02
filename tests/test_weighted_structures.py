
import pytest
from typing import List

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups,SingleEnsembleGroup, Sara2StructureList, Sara2SecondaryStructure, ComparisonStructures
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface, NupackSettings, MaterialParameter
from serena.weighted_structures import WeightedStructureResult, WeightedStructures,WeightedComparisonResult, WeightedGroupResult, MultipleGroupRawResults, IdealRangeSettings, PredictionResult

"""
def test_make_weighted_struct():
    nupack: NUPACK4Interface = NUPACK4Interface()
    sequence = 'GCCAUCGCAUGAGGAUAUGCUCCGGUUUCCGGAGCAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
    target = '........(((......(((.............))).....)))........................................'
    folded = '((((((.((((......((((((((...)))))))).....))))...(((.(((((.((....)))))))..))).)))))).'
    span = 5
    units = .5
    name = "09_eli"
    designID = 12345
    labname = "Tbox Round 1"
    folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/debug'
    ligand_oligo_energy:float = 10
    folded_energy_ligoligo: float = -29
    ligand_oligo_name:str = ''
    eterna_score:float = 100
    fold_change:float = 500
    number_of_clusters:int = 1000
    temp_C: int = 37
    rna_model: MaterialParameter = MaterialParameter.rna95_nupack4
    settings: NupackSettings = NupackSettings(material_param=rna_model, kcal_delta_span_from_mfe=span,
                                            temp_C=temp_C,Kcal_unit_increments=units,
                                            folded_2nd_state_structure=folded, folded_2nd_state_kcal=folded_energy_ligoligo,
                                            sequence=sequence)
    ensemble_groups: MultipleEnsembleGroups = MultipleEnsembleGroups()
    ensemble_groups = nupack.get_ensemble_groups(settings)
    single_group:SingleEnsembleGroup = ensemble_groups.groups[0]
    current_weight_struct: WeightedStructureResult = WeightedStructures.make_weighted_struct(single_group)
"""

def test_compair_weighted_structure():
    nupack: NUPACK4Interface = NUPACK4Interface()
    sequence = 'GCCAUCGCAUGAGGAUAUGCUCCGGUUUCCGGAGCAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
    target = '........(((......(((.............))).....)))........................................'
    folded = '((((((.((((......((((((((...)))))))).....))))...(((.(((((.((....)))))))..))).)))))).'
    span = 5
    units = .5
    name = "09_eli"
    designID = 12345
    labname = "Tbox Round 1"
    folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/debug'
    ligand_oligo_energy:float = 10
    folded_energy_ligoligo: float = -29
    ligand_oligo_name:str = ''
    eterna_score:float = 100
    fold_change:float = 500
    number_of_clusters:int = 1000
    temp_C: int = 37
    rna_model: MaterialParameter = MaterialParameter.rna95_nupack4
    settings: NupackSettings = NupackSettings(material_param=rna_model, kcal_delta_span_from_mfe=span,
                                            temp_C=temp_C,Kcal_unit_increments=units,
                                            folded_2nd_state_structure=folded, folded_2nd_state_kcal=folded_energy_ligoligo,
                                            sequence=sequence)
    ensemble_groups: MultipleEnsembleGroups = MultipleEnsembleGroups()
    ensemble_groups = nupack.get_ensemble_groups(settings)

    unbound_mfe_stuct:Sara2SecondaryStructure = Sara2SecondaryStructure(structure=ensemble_groups.non_switch_state_structure,
                                                                         freeEnergy=ensemble_groups.non_switch_state_mfe_kcal)
        
    bound_mfe_stuct:Sara2SecondaryStructure = Sara2SecondaryStructure(structure=ensemble_groups.switched_state_structure,
                                                                         freeEnergy=ensemble_groups.switched_state_mfe_kcal)

    single_group:SingleEnsembleGroup = ensemble_groups.groups[0]
    weighted : WeightedStructures = WeightedStructures()
    current_weight_struct: WeightedStructureResult = weighted.make_weighted_struct(single_group)
    nuc_count: int = current_weight_struct.weighted_struct.nuc_count
    current_comp_struct: WeightedComparisonResult = weighted.compair_weighted_structure(unbound_mfe_struct=unbound_mfe_stuct,
                                                                                    bound_mfe_struct=bound_mfe_stuct,
                                                                                    weighted_result=current_weight_struct,
                                                                                    nuc_count=nuc_count)
    assert nuc_count == 84
    assert current_comp_struct.num_both == 57







"""
Local Minima Structure Variation Data
Creation Date=2023-04-09 23:04:53.341105
---------------------------------------
***DESIGN INFO***
Design Name = Sara mod 1 of 6423414 #SaraFilterEverything-Good
DesignID = 6458932
Lab Name = Same State NG 1
Sequence = GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA
Eterna_Score = 100.0
FoldChange = 27.05
2nd State Target Structure = ........(((......(((.............))).....)))........................................
2nd State Folded Structure = ((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).
2nd State Folded Oligo Energy = -27.4
Energy Span from MFE = 7
Energy span units = .5
---------------------------------------
***RAW DATA***
Kcal,LMSV_U_mfe,LMSV_U_rel,LMSV_US_target,LMSV_US_folded
-31.679424285888672,1.0,1.0,48.0,22.0
-31.179424285888672,1.0,1.0,48.0,22.0
-30.679424285888672,2.0,0.0,51.0,23.0
-30.179424285888672,4.875,5.625,48.0,22.625
-29.679424285888672,4.571428571428571,4.857142857142858,47.0,22.71428571428571
-29.179424285888672,7.6875,8.0625,46.4375,17.625
-28.679424285888672,7.628571428571428,12.485714285714286,46.1142857142857,19.942857142857143
-28.179424285888672,8.92537313432836,19.95522388059703,46.67164179104479,20.014925373134336
-27.679424285888672,10.55882352941176,20.382352941176485,45.87254901960785,18.9607843137255
-27.179424285888672,10.880434782608695,19.266304347826082,45.88586956521739,19.30978260869564
-26.679424285888672,11.376811594202893,21.602898550724642,45.6608695652174,20.34202898550724
-26.179424285888672,12.073883161512022,22.18213058419245,45.687285223367695,20.096219931271488
-25.679424285888672,12.798024149286492,19.895718990120745,45.54774972557628,19.941822173435792
-25.179424285888672,13.47939262472885,18.55097613882864,45.09110629067244,19.798264642082437
---------------------------------------
***METRICS***
Polymorphicity Level (2kcal to end of sample) = 10.592338900810564
LMV_US_folded at folded with lignad/oligo energy = 19.30978260869564
---------------------------------------
EOF
"""