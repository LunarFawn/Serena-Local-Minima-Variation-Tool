
from scripts.run_switch_analysis import OriginalSwitchAnalysis, PredictionReponse, SwitchPrediction


def test_original_analysis():
    print("using debug")
    sequence = 'AUGGAUAUCACAGGAUAUGCUAUGGAUACCAUAGCAGAAGGGUGAUCCAUCCUAUGUAAAUACAUGAGGAUCACCCAUGUGUCC'
    target = '........(((......(((.............))).....)))........................................'
    folded = '....(((((....)))))(((((((...)))))))....(((((((((....(((((....))))).)))))))))........'
    folded_energy_ligoligo: float =  -30.6
    span = 7
    units = 1
    name = "09_eli"
    designID = 12345
    labname = "Tbox Round 1"
    folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/debug'
    ligand_oligo_energy:float = 10

    ligand_oligo_name:str = ''
    eterna_score:float = 100
    fold_change:float = 500
    number_of_clusters:int = 1000

    switch:OriginalSwitchAnalysis = OriginalSwitchAnalysis()

    analysis:PredictionReponse = switch.test_LMV(sequence=sequence,
                                                folded=folded,
                                                folded_energy_ligoligo=folded_energy_ligoligo,
                                                span=span,
                                                units=units,
                                                manual=False)

    print()
    print(analysis.message)
    print(f'Predicted Foldchage is {analysis.foldchange}')
    print(f'Predicted Switchyness is {analysis.prediction.name}')
    assert analysis.foldchange == -36.56982857142857
    assert analysis.prediction == SwitchPrediction.BAD