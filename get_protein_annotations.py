""" This file contains code for retreiving the experimentally determined
GO annotation for a list of proteins. """
from collections import defaultdict
from typing import Iterable, Optional, Generator
from Bio import SwissProt
from utils import (
    parse_go,
    has_go_annotation,
    NAMESPACES,
    EXPERIMENTAL_EVIDENCE_CODES,
)


def filter_swissprot_by_taxon(
    swissprot_parser: Generator, taxon_id: str
) -> Generator[SwissProt.Record, None, None]:
    """Generator function that filters a given instance of Bio.SwissProt.parse
    based on the given taxon_id.

    Yields Bio.SwissProt.Record instances
    """
    taxon_id = str(taxon_id)

    for sprot_record in swissprot_parser:
        if taxon_id in sprot_record.taxonomy_id:
            yield sprot_record


def filter_proteins(
    swissprot_filepath: str,
    taxon_id: Optional[str] = None,
    namespace_filter: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
) -> Generator[tuple, None, None]:
    """Generator function that yields two-element tuples of (1) protein_ids
    and (2) lists of go_terms that meet the given parameters (taxon,
    namespace_filter and allowed_evidence_codes)."""
    if allowed_evidence_codes is None:
        allowed_evidence_codes = EXPERIMENTAL_EVIDENCE_CODES

    with open(swissprot_filepath, "r") as swissprot_handle:
        swissprot_parser = SwissProt.parse(swissprot_handle)

        if taxon_id is not None:
            swissprot_records = filter_swissprot_by_taxon(
                swissprot_parser=swissprot_parser, taxon_id=taxon_id
            )
        else:
            swissprot_records = swissprot_parser

        for swissprot_record in swissprot_records:
            go_terms = parse_go(
                swissprot_record,
                namespace_filter=namespace_filter,
                evidence_code_filter=allowed_evidence_codes,
            )
            yield (swissprot_record.entry_name, go_terms)


def get_annotated_proteins(
    swissprot_filepath: str,
    taxon_id: Optional[str] = None,
    namespace_filter: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
) -> Generator[tuple, None, None]:
    """Generator function that yields tuples of protein IDs and lists of associated
    GO terms from the given SwissProt filepath (str)
    matching the given filtering parameters
    (taxon_id, namespace_filter, allowed_evidence_codes)
    """
    for protein, go_terms in filter_proteins(
        swissprot_filepath=swissprot_filepath,
        taxon_id=taxon_id,
        namespace_filter=namespace_filter,
        allowed_evidence_codes=allowed_evidence_codes,
    ):

        if len(go_terms) == 0:
            continue

        yield (protein, go_terms)


# TODO: Refactor this to use DRY code from above!
def get_proteins_under_annotation_threshold(
    sp_filepath: str,
    taxon_id: str,
    max_annotations: int = 3,
    namespace_filter: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
):
    """Generator function that yields two-element tuples of (1) protein_ids
    and (2) lists of go_terms that meet the given parameters (taxon,
    max_annotations, namespace_filter and allowed_evidence_codes)."""
    if allowed_evidence_codes is None:
        allowed_evidence_codes = EXPERIMENTAL_EVIDENCE_CODES

    sp_handle = open(sp_filepath, "r")

    for sp_record in SwissProt.parse(sp_handle):

        if str(taxon_id) not in sp_record.taxonomy_id:
            continue

        go_terms = parse_go(
            sp_record,
            namespace_filter=namespace_filter,
            evidence_code_filter=allowed_evidence_codes,
        )

        if len(go_terms) > max_annotations:
            continue

        yield (sp_record.entry_name, go_terms)

    sp_handle.close()


def get_annotated_protein_counts_by_namespace(
    sp_handle,
    taxons: Iterable[str],
    namespaces: Optional[Iterable] = NAMESPACES,
    evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
) -> dict:
    """For each GO namespace, number of proteins annotated with allowed
    terms. IMPORTANT TO NOTE: This counts proteins with any annotation for the
    given namespace, NOT the actual number of annotations.

    Uses experimental evidence codes, by default.

    Returns a dict with complex keys of (taxonomy, ontology namespace) and
    values of protein counts with annotation for the key pair.
    """
    taxons_set = {str(taxon) for taxon in taxons}
    counts_by_namespace = defaultdict(int)
    sp_handle.seek(0)
    sprot = SwissProt.parse(sp_handle)
    # Filter the sprot data to only the taxonomies of interest:
    sprot_records = [
        rec
        for rec in sprot
        if has_go_annotation(rec) and bool(taxons_set & set(rec.taxonomy_id))
    ]

    for sp_record in sprot_records:
        # sp_record.taxonomy_id is a list:
        tax_intersect = taxons_set.intersection(set(sp_record.taxonomy_id))
        # Is this valid? Can we assume there's only one taxon here?
        _taxon = list(tax_intersect)[0]

        go_list = parse_go(
            sp_record,
            namespace_filter=namespaces,
            evidence_code_filter=evidence_codes,
        )

        for namespace in namespaces:
            counts_by_namespace[(_taxon, namespace)] += int(
                any([term.ontology in namespace for term in go_list])
            )

    return counts_by_namespace


def count_annotated_proteins_from_files(
    sprot_file_list: Iterable,
    taxons: Iterable[str],
    allowed_namespaces: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
) -> dict:
    """Count the number of proteins with experimental annotations for the
    given taxonomy, namespaces and evidence codes in each of the given swissprot .dat files
    in sprot_file_list

    Returns a dictionary mapping filenames (keys) to protein counts (values)
    for the given taxonomy.
    """
    taxons = [str(taxon_id) for taxon_id in taxons]
    sp_count_by_namespace = []

    for filename in sprot_file_list:
        with open(filename, "r") as sp_handle:

            annotated_protein_count = get_annotated_protein_counts_by_namespace(
                sp_handle,
                taxons,
                namespaces=allowed_namespaces,
                evidence_codes=allowed_evidence_codes,
            )

            sp_count_by_namespace.append(
                {"filename": filename, "counts": annotated_protein_count}
            )

    return sp_count_by_namespace


def print_annotation_counts_table(data: dict, taxon: str):
    """ rudimentary pretty printing of annotation counts per ontology namespace. """

    print(data)
    print(taxon)

    # This is for getting column widths correct:
    col1_width = max([len(res.get("filename")) for res in data])

    # Get the ontologies as the keys of the first nested dict:
    namespaces_set = {k[1] for k in [list(s.get("counts").keys()) for s in data][0]}
    print(namespaces_set)

    preheader = f"| COUNTS OF PROTEINS WITH ANNOTATION {' ' * (col1_width+10)}|"
    table_header = "| FILE{} | TAXON | {}    | {}    | {}    | GROWTH |"
    table_row_sep = (
        f"{'-' * col1_width}------------------------------------------------"
    )
    table_row_pattern = "| {} | {}  | {:6,} | {:6,} | {:6,} | {} |"

    print(table_row_sep)
    print(preheader)
    print(table_row_sep)
    print(table_header.format(" " * (col1_width - 4), *namespaces_set))
    print(table_row_sep)
    previous_row_sum = None

    counts = [row.get("counts") for row in data]

    print(counts)

    for row in counts:
        # row should be a dict:
        for key, row_count in row.items():

            taxon, namespace = key
            print(taxon, namespace, row_count)
            row_sum = sum(row.values())

            if previous_row_sum is None:
                total_growth = " " * 6
            else:
                total_growth = f"{(row_sum-previous_row_sum):+6,}"

            _row = table_row_pattern.format(key, taxon, *row.values(), total_growth)

            print(_row)
            print(table_row_sep)
            previous_row_sum = row_sum


def get_annotations(
    protein_list_filepath: str,
    sprot_filepath: str,
    namespace_filter: Iterable,
    evidence_code_filter: Iterable,
):
    """This is intended to identify CAFA benchmark proteins by taking a list of
    protein IDs and identifying any GO annotations in a given swissprot file for
    the proteins in that list.

    Takes
    (1) filepath expects it to be a list of protein IDs, one per line
    (2) filepath for a Swissprot .dat file
    (3) GO namespaces to restrict search to
    (4) evidence codes to restrict search to

    outputs annotations (where they exist) for the given proteins
    """
    with open(protein_list_filepath, "r") as proteins_handle, open(
        sprot_filepath, "r"
    ) as sprot_handle:

        proteins_of_interest = [line.rstrip() for line in proteins_handle]
        sprot = SwissProt.parse(sprot_handle)

        for record in sprot:
            if record.entry_name in proteins_of_interest:
                go_terms = parse_go(
                    record,
                    namespace_filter=namespace_filter,
                    evidence_code_filter=evidence_code_filter,
                )

                for association in go_terms:
                    yield association


def print_annotations(
    protein_list_filepath: str,
    sprot_filepath: str,
    namespace_filter: Iterable,
    evidence_code_filter: Iterable,
):
    """This is a thin wrapper around get_annotations() that simple prints
    the data from get_annotations() to stdout.

    It can be used to generate groundtruth data for CAFA evaluation phases.

    TODO: However we still need to consider the "leaf-only" aspect of the benchmark.
    How do we filter this dataset to "leaf-only"?
    """

    annotations = get_annotations(
        protein_list_filepath,
        sprot_filepath,
        namespace_filter=namespace_filter,
        evidence_code_filter=evidence_code_filter,
    )
    for annotation in annotations:
        print("\t".join((annotation.protein, annotation.go_id)))

    annotations.close()


if __name__ == "__main__":
    PROTEINS_FILEPATH = "./data_for_blast/v2/9606_proteins_no_annotation.txt"
    SPROT_FILEPATH = (
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot.dat"
    )
    NAMESPACES = ("CCO",)
    EVIDENCE_CODES = EXPERIMENTAL_EVIDENCE_CODES
    TAXON_ID = 9606

    # THIS IS HOW ONE MIGHT GENERATE A LIST OF UNANNOTATED CANDIDATE PROTEINS:
    """
    proteins_without_annotation = get_proteins_under_annotation_threshold(
        SPROT_FILEPATH, TAXON_ID, max_annotations=0, namespace_filter=NAMESPACES
    )

    for protein, go_terms in test:
        print(protein)
    """
    # Alternatively:
    with open(PROTEINS_FILEPATH, "w") as write_handle:
        proteins_without_annotation = get_proteins_under_annotation_threshold(
            SPROT_FILEPATH, TAXON_ID, max_annotations=0, namespace_filter=NAMESPACES
        )
        # Here, the results are two-element tuples. The second element is
        # go_terms and will be an empty list that we can ignore:
        for protein, _ in proteins_without_annotation:
            write_handle.write(f"{protein}\n")

    # THIS IS HOW ONE WOULD COLLECT THE EXPERIMENTAL GROWTH FOR A LIST OF PROTEINS:
    # print_annotations(proteins_filepath, sprot_filepath, namespaces, evidence_codes)
