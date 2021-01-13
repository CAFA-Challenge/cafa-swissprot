import pytest
from filter_sprot_species import generate_protein_ids_mapping, get_no_exp_proteins


def test_unique_ids():
    ''' Test that all of the proteins yielded by get_no_exp_proteins() have
    unique generated IDs'''

    sprot_filepath = "./tests/test_data/uniprot/2019_01/uniprot_sprot_truncated.dat"
    taxonomies = (9606, 4577)
    ontologies = ("CCO", "BPO", "MFO")
    allowed_evidence_codes = (
        "IEA",
        "NR",
        "ND",
        "IC",
        "NAS",
        "TAS",
        "ISS",
        "ISO",
        "ISA",
        "ISM",
        "IGC",
        "IBA",
        "IBD",
        "IKR",
        "IRD",
        "RCA",
    )

    with open(sprot_filepath, "r") as sprot_handle:
        for taxon in taxonomies:
            proteins = get_no_exp_proteins(
                sprot_file_handle=sprot_handle,
                taxonomy=taxon,
                namespaces=ontologies,
                allowed_evidence_codes=allowed_evidence_codes,
            )

            cafa_ids = []
            sprot_ids = []

            for protein, namespace_map in proteins:
                cafa_ids.append(protein.id)
                sprot_ids.append(protein.description)

            assert len(cafa_ids) == len(sprot_ids)

            # If there are duplicate IDs in either list, that will be revealed
            # by comparing the lists to sets generated from the lists themselves:
            cafa_ids_unique = set(cafa_ids)
            sprot_ids_unique = set(sprot_ids)

            assert len(cafa_ids) == len(cafa_ids_unique)
            assert len(sprot_ids) == len(sprot_ids_unique)


if __name__ == "__main__":
    test_unique_ids()