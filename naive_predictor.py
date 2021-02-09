""" Code for generating a naive CAFA predictor. """
import pandas as pd
from matrix_utils import get_binary_annotation_matrix

def generate_naive_predictor(
        annotation_filepath: str, obo_filepath: str, namespace_verbose: str
) -> pd.Series:
    """ Returns a pandas.Series mapping GO term IDs to "predicted" probabilities """

    annotation_propagated_matrix = get_binary_annotation_matrix(
        annotation_filepath=annotation_filepath,
        obo_filepath=obo_filepath,
        namespace_verbose=namespace_verbose,
    )
    # Construct a Series of the probabilities from the binary Dataframe:
    naive_probabilities = (
            annotation_propagated_matrix.sum(axis=0) / annotation_propagated_matrix.shape[0]
    )
    return naive_probabilities


def main():
    """
    Generate a naive predictor and validate it's values against know Matlab values
    For the same inputs.

    Example usage and testing against Matlab-generated data. """
    from math import isclose
    annotation_filepath = "./input_annotation_for_blast_predictor.txt"
    obo_filepath = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    namespace = "cellular_component"

    naive_predictor = generate_naive_predictor(
        annotation_filepath=annotation_filepath,
        obo_filepath=obo_filepath,
        namespace_verbose=namespace,
    )


    print(naive_predictor.head(n=20))

    print("RUNNING BASIC ASSERTION TESTS")

    absolute_tolerance = 0.0001
    # 400 Sample values taken from a Matlab naive predictor generated from the
    # same annotation and obo files:
    test_data_from_matlab = (
        ("GO:0000015", 0.00000),
        ("GO:0000109", 0.00053),
        ("GO:0000110", 0.00000),
        ("GO:0000111", 0.00000),
        ("GO:0000112", 0.00000),
        ("GO:0000113", 0.00000),
        ("GO:0000118", 0.00133),
        ("GO:0000120", 0.00013),
        ("GO:0000123", 0.00373),
        ("GO:0000124", 0.00027),
        ("GO:0000125", 0.00013),
        ("GO:0000126", 0.00013),
        ("GO:0000127", 0.00067),
        ("GO:0000131", 0.00000),
        ("GO:0000133", 0.00000),
        ("GO:0000136", 0.00000),
        ("GO:0000137", 0.00000),
        ("GO:0000138", 0.00013),
        ("GO:0000139", 0.00800),
        ("GO:0000142", 0.00000),
        ("GO:0000144", 0.00000),
        ("GO:0000145", 0.00027),
        ("GO:0000148", 0.00000),
        ("GO:0000151", 0.01279),
        ("GO:0000152", 0.00200),
        ("GO:0000153", 0.00027),
        ("GO:0000159", 0.00067),
        ("GO:0000164", 0.00013),
        ("GO:0000172", 0.00040),
        ("GO:0000176", 0.00013),
        ("GO:0000177", 0.00000),
        ("GO:0000178", 0.00067),
        ("GO:0000214", 0.00000),
        ("GO:0000220", 0.00000),
        ("GO:0000221", 0.00000),
        ("GO:0000222", 0.00000),
        ("GO:0000223", 0.00000),
        ("GO:0000228", 0.01306),
        ("GO:0000229", 0.00013),
        ("GO:0000235", 0.00000),
        ("GO:0000242", 0.00013),
        ("GO:0000243", 0.00000),
        ("GO:0000262", 0.00000),
        ("GO:0000274", 0.00000),
        ("GO:0000275", 0.00013),
        ("GO:0000276", 0.00000),
        ("GO:0000306", 0.00013),
        ("GO:0000307", 0.00093),
        ("GO:0000308", 0.00013),
        ("GO:0000311", 0.00000),
        ("GO:0000312", 0.00000),
        ("GO:0000313", 0.00853),
        ("GO:0000314", 0.00200),
        ("GO:0000315", 0.00640),
        ("GO:0000322", 0.00000),
        ("GO:0000323", 0.01133),
        ("GO:0000324", 0.00000),
        ("GO:0000325", 0.00000),
        ("GO:0000326", 0.00000),
        ("GO:0000327", 0.00000),
        ("GO:0000328", 0.00000),
        ("GO:0000329", 0.00000),
        ("GO:0000330", 0.00000),
        ("GO:0000331", 0.00000),
        ("GO:0000333", 0.00000),
        ("GO:0000343", 0.00000),
        ("GO:0000344", 0.00000),
        ("GO:0000345", 0.00000),
        ("GO:0000346", 0.00053),
        ("GO:0000347", 0.00027),
        ("GO:0000399", 0.00000),
        ("GO:0000407", 0.00067),
        ("GO:0000408", 0.00013),
        ("GO:0000417", 0.00000),
        ("GO:0000418", 0.00000),
        ("GO:0000419", 0.00000),
        ("GO:0000421", 0.00053),
        ("GO:0000427", 0.00000),
        ("GO:0000428", 0.00493),
        ("GO:0000438", 0.00000),
        ("GO:0000439", 0.00013),
        ("GO:0000440", 0.00000),
        ("GO:0000444", 0.00027),
        ("GO:0000445", 0.00027),
        ("GO:0000446", 0.00000),
        ("GO:0000500", 0.00000),
        ("GO:0000502", 0.00426),
        ("GO:0000506", 0.00093),
        ("GO:0000775", 0.00400),
        ("GO:0000776", 0.00307),
        ("GO:0000777", 0.00080),
        ("GO:0000778", 0.00000),
        ("GO:0000779", 0.00093),
        ("GO:0000780", 0.00013),
        ("GO:0000781", 0.00067),
        ("GO:0000782", 0.00013),
        ("GO:0000783", 0.00013),
        ("GO:0000784", 0.00040),
        ("GO:0000785", 0.01359),
        ("GO:0000786", 0.00160),
        ("GO:0000787", 0.00000),
        ("GO:0000788", 0.00080),
        ("GO:0000789", 0.00013),
        ("GO:0000790", 0.00999),
        ("GO:0000791", 0.00080),
        ("GO:0000792", 0.00200),
        ("GO:0000793", 0.00253),
        ("GO:0000794", 0.00040),
        ("GO:0000795", 0.00013),
        ("GO:0000796", 0.00040),
        ("GO:0000797", 0.00000),
        ("GO:0000798", 0.00000),
        ("GO:0000799", 0.00013),
        ("GO:0000800", 0.00000),
        ("GO:0000801", 0.00000),
        ("GO:0000802", 0.00000),
        ("GO:0000803", 0.00040),
        ("GO:0000804", 0.00000),
        ("GO:0000805", 0.00040),
        ("GO:0000806", 0.00000),
        ("GO:0000807", 0.00000),
        ("GO:0000808", 0.00027),
        ("GO:0000809", 0.00000),
        ("GO:0000811", 0.00000),
        ("GO:0000812", 0.00027),
        ("GO:0000813", 0.00027),
        ("GO:0000814", 0.00000),
        ("GO:0000815", 0.00027),
        ("GO:0000817", 0.00000),
        ("GO:0000818", 0.00000),
        ("GO:0000835", 0.00013),
        ("GO:0000836", 0.00013),
        ("GO:0000837", 0.00000),
        ("GO:0000838", 0.00000),
        ("GO:0000839", 0.00013),
        ("GO:0000922", 0.00253),
        ("GO:0000923", 0.00000),
        ("GO:0000924", 0.00000),
        ("GO:0000927", 0.00000),
        ("GO:0000928", 0.00000),
        ("GO:0000930", 0.00040),
        ("GO:0000931", 0.00040),
        ("GO:0000932", 0.00200),
        ("GO:0000933", 0.00000),
        ("GO:0000934", 0.00000),
        ("GO:0000935", 0.00000),
        ("GO:0000936", 0.00000),
        ("GO:0000937", 0.00000),
        ("GO:0000938", 0.00000),
        ("GO:0000939", 0.00000),
        ("GO:0000940", 0.00040),
        ("GO:0000941", 0.00000),
        ("GO:0000942", 0.00000),
        ("GO:0000943", 0.00000),
        ("GO:0000974", 0.00000),
        ("GO:0001114", 0.00000),
        ("GO:0001400", 0.00000),
        ("GO:0001401", 0.00000),
        ("GO:0001405", 0.00000),
        ("GO:0001411", 0.00000),
        ("GO:0001518", 0.00080),
        ("GO:0001520", 0.00000),
        ("GO:0001527", 0.00053),
        ("GO:0001533", 0.00187),
        ("GO:0001534", 0.00000),
        ("GO:0001535", 0.00000),
        ("GO:0001536", 0.00000),
        ("GO:0001650", 0.00693),
        ("GO:0001651", 0.00000),
        ("GO:0001652", 0.00040),
        ("GO:0001669", 0.00133),
        ("GO:0001673", 0.00013),
        ("GO:0001674", 0.00000),
        ("GO:0001725", 0.00107),
        ("GO:0001726", 0.00213),
        ("GO:0001739", 0.00000),
        ("GO:0001740", 0.00040),
        ("GO:0001741", 0.00000),
        ("GO:0001750", 0.00093),
        ("GO:0001772", 0.00093),
        ("GO:0001891", 0.00040),
        ("GO:0001917", 0.00120),
        ("GO:0001931", 0.00040),
        ("GO:0001939", 0.00000),
        ("GO:0001940", 0.00000),
        ("GO:0002079", 0.00000),
        ("GO:0002080", 0.00013),
        ("GO:0002081", 0.00000),
        ("GO:0002095", 0.00000),
        ("GO:0002096", 0.00000),
        ("GO:0002102", 0.00040),
        ("GO:0002111", 0.00000),
        ("GO:0002116", 0.00013),
        ("GO:0002133", 0.00000),
        ("GO:0002139", 0.00000),
        ("GO:0002140", 0.00000),
        ("GO:0002141", 0.00000),
        ("GO:0002142", 0.00000),
        ("GO:0002144", 0.00000),
        ("GO:0002167", 0.00000),
        ("GO:0002169", 0.00000),
        ("GO:0002177", 0.00000),
        ("GO:0002178", 0.00080),
        ("GO:0002179", 0.00000),
        ("GO:0002180", 0.00000),
        ("GO:0002185", 0.00000),
        ("GO:0002186", 0.00000),
        ("GO:0002187", 0.00000),
        ("GO:0002189", 0.00000),
        ("GO:0002193", 0.00013),
        ("GO:0002197", 0.00000),
        ("GO:0002199", 0.00000),
        ("GO:0002929", 0.00000),
        ("GO:0002944", 0.00027),
        ("GO:0002945", 0.00013),
        ("GO:0002947", 0.00000),
        ("GO:0005575", 1.00000),
        ("GO:0005576", 0.06423),
        ("GO:0005577", 0.00013),
        ("GO:0005578", 0.00333),
        ("GO:0005579", 0.00080),
        ("GO:0005581", 0.00173),
        ("GO:0005582", 0.00000),
        ("GO:0005583", 0.00107),
        ("GO:0005584", 0.00027),
        ("GO:0005585", 0.00013),
        ("GO:0005586", 0.00013),
        ("GO:0005587", 0.00027),
        ("GO:0005588", 0.00027),
        ("GO:0005589", 0.00000),
        ("GO:0005590", 0.00000),
        ("GO:0005591", 0.00000),
        ("GO:0005592", 0.00013),
        ("GO:0005593", 0.00040),
        ("GO:0005594", 0.00040),
        ("GO:0005595", 0.00000),
        ("GO:0005596", 0.00000),
        ("GO:0005597", 0.00000),
        ("GO:0005598", 0.00000),
        ("GO:0005599", 0.00000),
        ("GO:0005600", 0.00000),
        ("GO:0005601", 0.00000),
        ("GO:0005602", 0.00000),
        ("GO:0005604", 0.00173),
        ("GO:0005605", 0.00027),
        ("GO:0005606", 0.00000),
        ("GO:0005607", 0.00000),
        ("GO:0005608", 0.00013),
        ("GO:0005609", 0.00000),
        ("GO:0005610", 0.00000),
        ("GO:0005611", 0.00000),
        ("GO:0005612", 0.00000),
        ("GO:0005614", 0.00000),
        ("GO:0005615", 0.04611),
        ("GO:0005616", 0.00000),
        ("GO:0005618", 0.00000),
        ("GO:0005619", 0.00000),
        ("GO:0005621", 0.00000),
        ("GO:0005622", 0.80677),
        ("GO:0005623", 0.92417),
        ("GO:0005628", 0.00000),
        ("GO:0005630", 0.00000),
        ("GO:0005631", 0.00000),
        ("GO:0005632", 0.00000),
        ("GO:0005633", 0.00000),
        ("GO:0005634", 0.35181),
        ("GO:0005635", 0.01333),
        ("GO:0005637", 0.00080),
        ("GO:0005638", 0.00000),
        ("GO:0005639", 0.00027),
        ("GO:0005640", 0.00027),
        ("GO:0005641", 0.00013),
        ("GO:0005642", 0.00013),
        ("GO:0005643", 0.00187),
        ("GO:0005652", 0.00000),
        ("GO:0005654", 0.22881),
        ("GO:0005655", 0.00000),
        ("GO:0005656", 0.00000),
        ("GO:0005657", 0.00293),
        ("GO:0005658", 0.00013),
        ("GO:0005662", 0.00027),
        ("GO:0005663", 0.00053),
        ("GO:0005664", 0.00013),
        ("GO:0005665", 0.00133),
        ("GO:0005666", 0.00093),
        ("GO:0005667", 0.00733),
        ("GO:0005668", 0.00000),
        ("GO:0005669", 0.00093),
        ("GO:0005671", 0.00093),
        ("GO:0005672", 0.00000),
        ("GO:0005673", 0.00000),
        ("GO:0005674", 0.00000),
        ("GO:0005675", 0.00013),
        ("GO:0005677", 0.00013),
        ("GO:0005680", 0.00187),
        ("GO:0005681", 0.00560),
        ("GO:0005682", 0.00000),
        ("GO:0005683", 0.00000),
        ("GO:0005684", 0.00160),
        ("GO:0005685", 0.00027),
        ("GO:0005686", 0.00000),
        ("GO:0005687", 0.00000),
        ("GO:0005688", 0.00000),
        ("GO:0005689", 0.00120),
        ("GO:0005690", 0.00000),
        ("GO:0005691", 0.00000),
        ("GO:0005692", 0.00000),
        ("GO:0005693", 0.00000),
        ("GO:0005694", 0.02652),
        ("GO:0005697", 0.00053),
        ("GO:0005700", 0.00000),
        ("GO:0005701", 0.00000),
        ("GO:0005702", 0.00000),
        ("GO:0005703", 0.00000),
        ("GO:0005704", 0.00000),
        ("GO:0005705", 0.00000),
        ("GO:0005706", 0.00000),
        ("GO:0005712", 0.00000),
        ("GO:0005713", 0.00000),
        ("GO:0005714", 0.00000),
        ("GO:0005715", 0.00000),
        ("GO:0005719", 0.00080),
        ("GO:0005720", 0.00120),
        ("GO:0005721", 0.00040),
        ("GO:0005722", 0.00000),
        ("GO:0005723", 0.00000),
        ("GO:0005724", 0.00000),
        ("GO:0005725", 0.00000),
        ("GO:0005726", 0.00000),
        ("GO:0005727", 0.00000),
        ("GO:0005728", 0.00000),
        ("GO:0005729", 0.00000),
        ("GO:0005730", 0.04558),
        ("GO:0005731", 0.00000),
        ("GO:0005732", 0.00093),
        ("GO:0005736", 0.00000),
        ("GO:0005737", 0.53545),
        ("GO:0005739", 0.09422),
        ("GO:0005740", 0.02492),
        ("GO:0005741", 0.00426),
        ("GO:0005742", 0.00040),
        ("GO:0005743", 0.01639),
        ("GO:0005744", 0.00027),
        ("GO:0005745", 0.00000),
        ("GO:0005746", 0.00600),
        ("GO:0005747", 0.00520),
        ("GO:0005749", 0.00000),
        ("GO:0005750", 0.00027),
        ("GO:0005751", 0.00080),
        ("GO:0005753", 0.00173),
        ("GO:0005754", 0.00000),
        ("GO:0005756", 0.00000),
        ("GO:0005757", 0.00027),
        ("GO:0005758", 0.00293),
        ("GO:0005759", 0.01692),
        ("GO:0005760", 0.00013),
        ("GO:0005761", 0.00853),
        ("GO:0005762", 0.00640),
        ("GO:0005763", 0.00200),
        ("GO:0005764", 0.01133),
        ("GO:0005765", 0.00520),
        ("GO:0005766", 0.00053),
        ("GO:0005767", 0.00053),
        ("GO:0005768", 0.02092),
        ("GO:0005769", 0.00533),
        ("GO:0005770", 0.00466),
        ("GO:0005771", 0.00013),
        ("GO:0005773", 0.03172),
        ("GO:0005774", 0.01119),
        ("GO:0005775", 0.00013),
        ("GO:0005776", 0.00267),
        ("GO:0005777", 0.00800),
        ("GO:0005778", 0.00173),
        ("GO:0005779", 0.00093),
        ("GO:0005780", 0.00000),
        ("GO:0005782", 0.00067),
        ("GO:0005783", 0.06903),
        ("GO:0005784", 0.00000),
        ("GO:0005785", 0.00000),
        ("GO:0005786", 0.00040),
        ("GO:0005787", 0.00000),
        ("GO:0005788", 0.00253),
        ("GO:0005789", 0.02305),
        ("GO:0005790", 0.00000),
        ("GO:0005791", 0.00053),
        ("GO:0005793", 0.00320),
        ("GO:0005794", 0.04837),
        ("GO:0005795", 0.00160),
        ("GO:0005796", 0.00040),
        ("GO:0005797", 0.00027),
        ("GO:0005798", 0.00133),
        ("GO:0005801", 0.00093),
        ("GO:0005802", 0.00586),
        ("GO:0005811", 0.00400),
        ("GO:0005813", 0.02785),
        ("GO:0005814", 0.00453),
        ("GO:0005815", 0.03038),
        ("GO:0005816", 0.00000),
        ("GO:0005818", 0.00000),
        ("GO:0005819", 0.00613),
    )

    for term, matlab_probability in test_data_from_matlab[:300]:
        print(naive_predictor.loc[term], matlab_probability)

        assert isclose(
            naive_predictor.loc[term], matlab_probability, abs_tol=absolute_tolerance
        )


if __name__ == "__main__":
    main()
