"""
pytest for ensemble variation
"""

import pytest

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups, Sara2StructureList, Sara2SecondaryStructure, ComparisonStructures
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface, NupackSettings, MaterialParameter
from serena.weighted_structures import WeightedStructures



#first setup nupack for folding

nupack: NUPACK4Interface = NUPACK4Interface()
sequence = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
target = '........(((......(((.............))).....)))........................................'
folded = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
span = 7
units = 1
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


rna_model: MaterialParameter = MaterialParameter.rna95_nupack3





settings: NupackSettings = NupackSettings(material_param=rna_model, kcal_delta_span_from_mfe=span,
                                          temp_C=temp_C,Kcal_unit_increments=units,
                                          folded_2nd_state_structure=folded, folded_2nd_state_kcal=folded_energy_ligoligo,
                                          sequence=sequence)

ensemble_groups: MultipleEnsembleGroups = nupack.get_ensemble_groups(settings)

comparison_structures:ComparisonStructures = ComparisonStructures()

#add unbound mfe
unbound_mfe_structure: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence, 
                                                                    structure=ensemble_groups.groups[0].group.mfe_structure, 
                                                                    freeEnergy=ensemble_groups.groups[0].group.mfe_freeEnergy)
mfe_struct_name:str = 'unbound'
comparison_structures.add_structure(unbound_mfe_structure, mfe_struct_name)


#add bound folded structures
folded_structure: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence, 
                                                                    structure=folded, 
                                                                    freeEnergy=folded_energy_ligoligo)
folded_struct_name:str = 'bound'
comparison_structures.add_structure(folded_structure, folded_struct_name)

ensemble_variation: EnsembleVariation = EnsembleVariation()
local_minima_variations_unbound: EVResult = ensemble_variation.get_ensemble_variation(ensemble=ensemble_groups, comparison_structure=comparison_structures.get_structure_by_name(mfe_struct_name))
local_minima_variations_bound: EVResult = ensemble_variation.get_ensemble_variation(ensemble=ensemble_groups, comparison_structure=comparison_structures.get_structure_by_name(folded_struct_name))


#now get weighted structure data

weighted : WeightedStructures = WeightedStructures()



print (local_minima_variations_bound)

print (local_minima_variations_bound)
