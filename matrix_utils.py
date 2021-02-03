from goatools.obo_parser import GODag
import pandas as pd

def get_binary_annotation_matrix(
        annotation_filepath: str, obo_filepath: str, namespace_verbose: str
) -> pd.DataFrame:
    """ Returns a binary pandas.DataFrame mapping protein IDs to GO term IDs"""

    optional_attrs = ["relationship", "replaced_by", "consider"]
    dag = GODag(
        obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None
    )

    # this is canonical IDs, not alt-IDs:
    canonical_namespace_terms = {
        term_obj.id: term_obj.get_all_upper()
        for term_id, term_obj in dag.items()
        if term_obj.namespace == namespace_verbose
    }

    # this list CAN and SHOULD contain alt-IDs:
    valid_namespace_terms = sorted(
        list(
            {
                term_id
                for term_id, term_object in dag.items()
                if term_object.namespace == namespace_verbose
            }
        )
    )

    with open(annotation_filepath, "r") as annotation_handle:
        # loop all the protein ID => GO term ID pairs in the input file and collect
        # all unique protein IDs where the corresponding GO term is valid:
        proteins = {
            protein
            for line in annotation_handle
            for protein, term in [line.rstrip().split()]
            if term in valid_namespace_terms
        }

        # A binary DataFrame with proteins (rows) and GO terms (columns):
        annotation_propagated_matrix = pd.DataFrame(
            0, index=sorted(proteins), columns=sorted(list(canonical_namespace_terms))
        )

        # Loop the file a second time to populate the DataFrame:
        annotation_handle.seek(0)

        for line in annotation_handle:

            protein, go_term = line.rstrip().split()
            if go_term not in valid_namespace_terms:
                continue

            lookup = dag.query_term(go_term)
            annotation_propagated_matrix.loc[protein, lookup.id] = 1
            # Propagate the annotation:
            annotation_propagated_matrix.loc[protein, lookup.get_all_upper()] = 1

    return annotation_propagated_matrix
