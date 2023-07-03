"""
test for weighted structures
"""

import pytest
from typing import List

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups, Sara2StructureList, Sara2SecondaryStructure, ComparisonStructures
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface, NupackSettings, MaterialParameter
from serena.weighted_structures import WeightedStructures, WeightedGroupResult, MultipleGroupRawResults, IdealRangeSettings, PredictionResult



#first setup nupack for folding

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

#now get weighted structure data

weighted : WeightedStructures = WeightedStructures()

temperatures: List[int] = [36,37,38]

weighted_raw_result:MultipleGroupRawResults = weighted.get_multi_temp_full_raw_data(temperature_list=temperatures, ensemble_groups=ensemble_groups)

ideal_ranges:IdealRangeSettings = IdealRangeSettings(bound_kcal_span_plus=3,
                                                     bound_kcal_span_minus=3,
                                                     mfe_effect_range_plus=2)


prediciton_results: List[PredictionResult] = weighted.find_ideal_switch_range_ensemble(raw_results=weighted_raw_result, settings=ideal_ranges)

print(prediciton_results)