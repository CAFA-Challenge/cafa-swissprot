import pytest
from experimental_growth import count_experimental_proteins, _experimental_evidence_codes


def test_simp():
    expected = {
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_01/uniprot_sprot_truncated.dat": {
            "BPO": 8786,
            "MFO": 7539,
            "CCO": 9905,
        },
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_06/uniprot_sprot_truncated.dat": {
            "BPO": 8960,
            "MFO": 7533,
            "CCO": 9996,
        },
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_01/uniprot_sprot_truncated.dat": {
            "BPO": 9091,
            "MFO": 7638,
            "CCO": 10153,
        },
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot_truncated.dat": {
            "BPO": 9282,
            "MFO": 8175,
            "CCO": 10294,
        },
    }
    
    '''
    ----------------------------------------------------------------------------------------------------------------------
| /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_01/uniprot_sprot_truncated.dat |  7,539 |  8,786 |  9,905 |
----------------------------------------------------------------------------------------------------------------------
| /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_06/uniprot_sprot_truncated.dat |  7,533 |  8,960 |  9,996 |
----------------------------------------------------------------------------------------------------------------------
| /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_01/uniprot_sprot_truncated.dat |  7,638 |  9,091 | 10,153 |
----------------------------------------------------------------------------------------------------------------------
| /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot_truncated.dat |  8,175 |  9,282 | 10,294 |
----------------------------------------------------------------------------------------------------------------------

    
    '''



    filelist = expected.keys()

    # python experimental_growth.py -f ../repo2/scott_uniprot_list.txt --taxon 9606
    test = count_experimental_proteins(filelist, 9606, _experimental_evidence_codes)

    for key, count_dict in test.items():
        assert expected.get(key) == count_dict

