import pytest



from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, NupackSettings, MaterialParameter
from serena.utilities.ensemble_structures import  Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs
from serena.local_minima_variation import LocalMinimaVariation
from serena.utilities.ensemble_variation import EVResult

@pytest.fixture
def nupack_4_settings():
    return NupackSettings()

@pytest.fixture
def initialized_nupack_4_settings():
    return NupackSettings(material_param=MaterialParameter.rna95_nupack4,
                          temp_C=37,
                          kcal_span_from_mfe=5,
                          Kcal_unit_increments=1,
                          sequence='ACGUACAUGAC')

@pytest.fixture
def nupack_switch_states():
    unbound_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence='ACGUACAUGAC',
                                                                     structure='(((....).))')
    bound_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence='ACGUACAUGAC',
                                                                     structure='((.......))')
    
    switch_state:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=unbound_struct,
                                                                               switched_mfe_struct=bound_struct)
    return switch_state
    


def test_empty_nupack_4_settings(nupack_4_settings:NupackSettings):
    assert nupack_4_settings.material_param == MaterialParameter.NONE
    assert nupack_4_settings.temp_C == 0
    assert nupack_4_settings.kcal_span_from_mfe == 0
    assert nupack_4_settings.Kcal_unit_increments == 0
    assert nupack_4_settings.sequence == ''

def test_initialized_nupack_4_settings(initialized_nupack_4_settings:NupackSettings):
    assert initialized_nupack_4_settings.material_param == MaterialParameter.rna95_nupack4
    assert initialized_nupack_4_settings.temp_C == 37
    assert initialized_nupack_4_settings.kcal_span_from_mfe == 5
    assert initialized_nupack_4_settings.Kcal_unit_increments == 1
    assert initialized_nupack_4_settings.sequence == 'ACGUACAUGAC'
    

def test_get_subopt_groups(initialized_nupack_4_settings:NupackSettings):
    nupack4: NUPACK4Interface = NUPACK4Interface()
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    assert structs.num_structures == 10
    assert structs.sara_stuctures[0].structure == '...........'
    assert structs.sara_stuctures[0].freeEnergy == 0.0
    assert structs.sara_stuctures[0].stackEnergy == 0.0
    assert structs.sara_stuctures[0].nuc_count == 11
    assert structs.sara_stuctures[0].sequence == 'ACGUACAUGAC'    
    
    assert structs.sara_stuctures[1].structure == '...(....)..'
    assert structs.sara_stuctures[1].freeEnergy == 2.0991673469543457
    assert structs.sara_stuctures[1].stackEnergy == 0.0
    assert structs.sara_stuctures[1].nuc_count == 11
    assert structs.sara_stuctures[1].sequence == 'ACGUACAUGAC'
     
    assert structs.sara_stuctures[2].structure == '...(.....).'
    assert structs.sara_stuctures[2].freeEnergy == 1.338149070739746
    assert structs.sara_stuctures[2].stackEnergy == 0.0
    assert structs.sara_stuctures[2].nuc_count == 11
    assert structs.sara_stuctures[2].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[3].structure == '..(....)...'
    assert structs.sara_stuctures[3].freeEnergy == 2.453690528869629
    assert structs.sara_stuctures[3].stackEnergy == 0.0
    assert structs.sara_stuctures[3].nuc_count == 11
    assert structs.sara_stuctures[3].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[4].structure == '..(.......)'
    assert structs.sara_stuctures[4].freeEnergy == 3.604722738265991
    assert structs.sara_stuctures[4].stackEnergy == 0.0
    assert structs.sara_stuctures[4].nuc_count == 11
    assert structs.sara_stuctures[4].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[5].structure == '..((....).)'
    assert structs.sara_stuctures[5].freeEnergy == 4.704723358154297
    assert structs.sara_stuctures[5].stackEnergy == 0.0
    assert structs.sara_stuctures[5].nuc_count == 11
    assert structs.sara_stuctures[5].sequence == 'ACGUACAUGAC'

    assert structs.sara_stuctures[6].structure == '..((.....))'
    assert structs.sara_stuctures[6].freeEnergy == -0.0952768325805664
    assert structs.sara_stuctures[6].stackEnergy == 0.0
    assert structs.sara_stuctures[6].nuc_count == 11
    assert structs.sara_stuctures[6].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[7].structure == '.(......)..'
    assert structs.sara_stuctures[7].freeEnergy == 1.1241339445114136
    assert structs.sara_stuctures[7].stackEnergy == 0.0
    assert structs.sara_stuctures[7].nuc_count == 11
    assert structs.sara_stuctures[7].sequence == 'ACGUACAUGAC'    
    
    assert structs.sara_stuctures[8].structure == '.((....))..'
    assert structs.sara_stuctures[8].freeEnergy == 1.2241339683532715
    assert structs.sara_stuctures[8].stackEnergy == 0.0
    assert structs.sara_stuctures[8].nuc_count == 11
    assert structs.sara_stuctures[8].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[9].structure == '(......)...'
    assert structs.sara_stuctures[9].freeEnergy == 3.1283416748046875
    assert structs.sara_stuctures[9].stackEnergy == 0.0
    assert structs.sara_stuctures[9].nuc_count == 11
    assert structs.sara_stuctures[9].sequence == 'ACGUACAUGAC'

def test_load_nupack_subopt_as_ensemble(initialized_nupack_4_settings:NupackSettings, nupack_switch_states:EnsembleSwitchStateMFEStructs):
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                             settings=initialized_nupack_4_settings,
                                                                             switch_state=nupack_switch_states
                                                                             )
    assert ensemble.num_groups == 5

def test_get_lmv_nupack_ensemble(initialized_nupack_4_settings:NupackSettings, nupack_switch_states:EnsembleSwitchStateMFEStructs):
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                             settings=initialized_nupack_4_settings,
                                                                             switch_state=nupack_switch_states
                                                                             )
    lmv:LocalMinimaVariation = LocalMinimaVariation()
    groups_results: EVResult = lmv.get_multi_group_lmv(ensemble=ensemble,
                                                       reference_structure=nupack_switch_states.non_switch_mfe_struct)
    assert groups_results.ev_values[0].ev_normalized == 6.0
    assert groups_results.ev_values[1].ev_normalized == 5.333333333333333
    assert groups_results.ev_values[2].ev_normalized == 6.0
    assert groups_results.ev_values[3].ev_normalized == 4.0
    assert groups_results.ev_values[4].ev_normalized == 6.0
    