import pytest



from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, NupackSettings, MaterialParameter



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
    nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )