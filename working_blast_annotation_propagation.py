import os
import sys
from typing import Iterable
import contextlib
import numpy as np
import pandas as pd
from Bio import SwissProt
from goatools.obo_parser import GODag
from utils import parse_go, EXPERIMENTAL_EVIDENCE_CODES


def get_propagation_map(obo_filepath:str, namespace:str) -> pd.DataFrame:
    """ Generates and returns an ontology annotation propagation map

    returns a pandas.DataFrame
    """
    #obo_filepath = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    optional_attrs = ["relationship", "replaced_by", "consider"]

    dag = GODag(
        obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None
    )

    # namespace = "CCO"
    #namespace_long = "cellular_component"

    propagation_map = {
        term_object.item_id: [*list(term_object.get_all_upper()), term_object.item_id]
        for term_id, term_object in dag.items()
        if term_object.namespace == namespace
    }

    df_keys = sorted(list(propagation_map.keys()))
    propagation_map_df = pd.DataFrame(data=0, index=df_keys, columns=df_keys)

    for key, val in propagation_map.items():
        propagation_map_df.loc[key, val] = 1

    #propagation_map_df.to_pickle("cco_propagation_map.pkl")
    return propagation_map_df


class BlastPredictorGenerator:
    """ For generating CAFA Blast predictors """

    def __init__(
        self,
        blastp_results_filepath: str,
        swissprot_filepath: str,
        propagation_map_filepath: str,
        obo_filepath: str,
    ):
        self.prediction_df = None
        self.blastp_results_filepath = blastp_results_filepath
        self.swissprot_filepath = swissprot_filepath
        self.propagation_map_df = pd.read_pickle(propagation_map_filepath)

        optional_dag_attrs = ["relationship", "replaced_by", "consider"]
        self.dag = GODag(
            obo_filepath,
            optional_attrs=optional_dag_attrs,
            load_obsolete=False,
            prt=None,
        )

    def get_predictions(
        self, query_protein_ids: Iterable, ontology_namespaces: Iterable
    ) -> pd.DataFrame:

        """Generates a Blast-based predictor protein -> GO term predictor
        based on sequence similarity and known annotation of those similar proteins.

        The queried proteins MUST be used as the query to the pre-build blastp results:
        blastp -query ./9606_no_annotation.fa -db data_for_blast/v2/uniprot_sprot_CCO
        -num_threads 6 -outfmt "6 qseqid sseqid evalue length pident nident"
        -out data_for_blast/v2/cafa_targets_blast_results.txt

        """
        # self.query_protein_ids = query_protein_ids
        # self.ontology_namespaces = ontology_namespaces

        # initialize an all-zero pandas DataFrame with protein as rows and ontology terms as cols:
        self.prediction_df = pd.DataFrame(
            data=0,
            index=query_protein_ids,
            columns=list(self.propagation_map_df.index),
        )

        # this gives us a dictionary mapping query protein IDs to tuples
        # containing similiar protein IDs and the similiarity scores:
        blastp_results = get_blastp_lookup_dict(
            query_protein_ids, self.blastp_results_filepath
        )

        # lookup all the similiar protein IDs in SwissProt and store them as a dictionary
        # mapping protein IDs (keys) to lists of GO annotations (values):
        swissprot_lookup_proteins = [
            blast_tuple[1]
            for nested_list in blastp_results.values()
            for blast_tuple in nested_list
        ]
        swissprot_lookup_dict = get_swissprot_lookup_dict(
            swissprot_lookup_proteins,
            self.swissprot_filepath,
            ontology_namespaces=ontology_namespaces,
        )

        # loop over all the input proteins and the corresponding similiar (via blastp)
        # proteins, identifying the annotation for the similiar proteins. Propagate
        # the annotations and assign to the input proteins using the % sequence
        # similarity as the confidence value
        for target_protein in query_protein_ids:
            for blast_hit in blastp_results[target_protein]:
                # blast_hit will be an Iterable of blastp results with this form:
                # ('ZZEF1_HUMAN', 'CUL7_HUMAN', 41.221)
                (
                    blast_hit_protein_id,
                    similarity_percentage,
                ) = blast_hit[1:]
                probability = np.round(similarity_percentage / 100, decimals=5)

                for leaf_annotation in swissprot_lookup_dict[blast_hit_protein_id]:
                    self.propagate_prediction(
                        target_protein, leaf_annotation.go_id, probability
                    )

        return self.prediction_df

    def propagate_prediction(
        self, protein_id: str, base_go_id: str, probability: float
    ):
        """Updates the row in self.prediction_df for the given protein
        with the given probability based on the propagation of the given GO ID.
        """
        try:
            propagation_df_row = self.propagation_map_df.loc[base_go_id, :]
        except KeyError:
            try:
                # Need to query on the go_id to see if it is an alt id
                # temporarily overriding sys.stderr.write in order to squash
                # annoying stderr messages from goatools.GODag.query_term()
                # when a term is not found:
                stash_stderr_write = sys.stderr.write
                sys.stderr.write = lambda err: 0
                go_lookup = self.dag.query_term(base_go_id)

                # restore sys.stderr.write
                sys.stderr.write = stash_stderr_write

                if go_lookup is not None:
                    go_id = go_lookup.item_id
                    propagation_df_row = self.propagation_map_df.loc[go_id, :]
                else:
                    return
            except KeyError:
                # the GO term ID is NOT in the propagation dataframe!
                # likely means it is obsolete...
                return

        propagation_mask = propagation_df_row == 1
        # column names (GO IDs) that are ancestors of our leaf GO ID:
        propagation_mask_cols = propagation_mask.index[propagation_mask].tolist()
        # The cells for the protein and the GO IDs of interest:
        propagation_cells = self.prediction_df.loc[protein_id, propagation_mask_cols]
        # Boolean pandas.Sequence representing if the existing confidence
        # float values are < the confidence value at hand:
        update_mask = propagation_cells < probability
        # the col names (GO IDs) to update
        # (have an existing val < the confidence at hand):
        cols_to_update = update_mask.index[update_mask].tolist()
        # Finally update the dataframe:
        self.prediction_df.loc[protein_id, cols_to_update] = probability


def get_swissprot_lookup_dict(
    lookup_protein_ids: Iterable,
    swissprot_filepath: str,
    ontology_namespaces: Iterable = (),
) -> dict:
    """Generates and returns a dict mapping protein IDs (keys) to lists of
    GO annotations.
    """
    swissprot_lookup_dict = {}  # protein_id: [] for protein_id in lookup_protein_ids}

    for swissprot_record in SwissProt.parse(swissprot_filepath):
        if swissprot_record.entry_name in lookup_protein_ids:
            swissprot_lookup_dict[swissprot_record.entry_name] = parse_go(
                swissprot_record,
                namespace_filter=ontology_namespaces,
                evidence_code_filter=EXPERIMENTAL_EVIDENCE_CODES,
            )
    return swissprot_lookup_dict


def get_blastp_lookup_dict(protein_ids: Iterable, blastp_results_filepath: str) -> dict:
    """Takes a filepath for blastp output and parses the file contents into
    a dict with query protein IDs as keys and tuples as values where the tuple
    has the form (query_protein_id, similar_protein_id, similarity_percentage)
    """
    blastp_results_dict = {protein: [] for protein in protein_ids}

    with open(blastp_results_filepath, "r") as blastp_results_handle:
        for line in blastp_results_handle:
            query_protein_id, *_ = line.rstrip().split()
            if query_protein_id in protein_ids:
                line_split = line.rstrip().split()
                similar_protein_id = line_split[1]
                similarity = float(line_split[4])
                blastp_results_dict[query_protein_id].append(
                    (query_protein_id, similar_protein_id, similarity)
                )
    return blastp_results_dict


if __name__ == "__main__":

    from tests.test_data.blast_predictor.per_protein_matlab_blast_prediction_scores import (
        test_blast_scores,
    )

    PROPAGATION_MAP_FILEPATH = "./cco_propagation_map.pkl"
    OBO_FILEPATH = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    SPROT_FILEPATH = (
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot.dat"
    )
    ONTOLOGIES = ("CCO",)


    BLASTP_SIMILARITY_FILEPATH = "./data_for_blast/v2/blastp.out"

    blast_generator = BlastPredictorGenerator(
        blastp_results_filepath=BLASTP_SIMILARITY_FILEPATH,
        swissprot_filepath=SPROT_FILEPATH,
        propagation_map_filepath=PROPAGATION_MAP_FILEPATH,
        obo_filepath=OBO_FILEPATH,
    )

    blast_predictions = blast_generator.get_predictions(
        query_protein_ids=set(test_blast_scores.keys()), ontology_namespaces=ONTOLOGIES
    )

    print(blast_predictions.head())

    # ====================================
    # Testing
    # ====================================
    # Test that the value of every non-zero cell in the returned DataFrame matches
    # the corresponding data harvested from Matlab output:
    print("TESTING GENERATED DATAFRAME AGAINST KNOWN DATA FROM MATLAB")

    for protein in blast_predictions.index:
        print(f"   TESTING {protein}")
        test_counter = 0
        for go_term in blast_predictions.loc[protein].index:
            # floating point number for protein/GO pair:
            v = blast_predictions.loc[protein, go_term]
            if v == 0:
                continue

            # correspondng floating point value from the Matlab data:
            test_v = test_blast_scores.get(protein).get("terms").get(go_term)
            assert v == test_v
            test_counter += 1

        print(f"      {test_counter} TESTS PASSED")
