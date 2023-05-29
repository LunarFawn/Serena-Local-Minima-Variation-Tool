"""
pytest for ensemble variation
"""

import pytest

import serena.ensemble_variation as ev
from serena.ensemble_variation import EnsembleVariation, EVResult
import serena.structures as ser_structs
from serena.structures import MultipleEnsembleGroups
import serena.nupack4_sara2_extension as nupack_extension
from serena.nupack4_sara2_extension import NUPACK4Interface



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



ensemble_groups: MultipleEnsembleGroups = nupack.get_ensemble_groups(sequence, kcal_delta_span_from_mfe=span, Kcal_unit_increments=units,
                                                          folded_2nd_state_structure=folded, folded_2nd_state_kcal=folded_energy_ligoligo)

ensemble_variation: EnsembleVariation = EnsembleVariation()


local_minima_variations_unbound: EVResult = ensemble_variation.get_ensemble_variation(ensemble=ensemble_groups, state_source=1)
local_minima_variations_bound: EVResult = ensemble_variation.get_ensemble_variation(ensemble=ensemble_groups, state_source=2)
