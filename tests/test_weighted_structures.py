
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
