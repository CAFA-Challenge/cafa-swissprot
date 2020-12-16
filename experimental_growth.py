from collections import defaultdict
from typing import Iterable, Optional
from Bio import SwissProt as sp
from utils import (
    parse_go,
    has_go_annotation,
    NAMESPACES,
    EXPERIMENTAL_EVIDENCE_CODES,
)


def get_proteins_under_annotation_threshold(
    sp_filepath: str,
    taxon_id: str,
    max_annotations: int = 3,
    namespace_filter: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
):

    if allowed_evidence_codes is None:
        allowed_evidence_codes = EXPERIMENTAL_EVIDENCE_CODES

    sp_handle = open(sp_filepath, "r")

    for sp_record in sp.parse(sp_handle):

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
    taxons = set([str(taxon) for taxon in taxons])
    counts_by_namespace = defaultdict(int)
    sp_handle.seek(0)
    sprot = sp.parse(sp_handle)
    # Filter the sprot data to only the taxonomies of interest:
    sprot_records = [
        rec
        for rec in sprot
        if has_go_annotation(rec) and bool(taxons & set(rec.taxonomy_id))
    ]

    for sp_record in sprot_records:
        # sp_record.taxonomy_id is a list:
        tax_intersect = taxons.intersection(set(sp_record.taxonomy_id))
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
    sp_file_list: Iterable,
    taxons: Iterable[str],
    allowed_namespaces: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = EXPERIMENTAL_EVIDENCE_CODES,
) -> dict:
    """Count the number of proteins with experimental annotations for the
    given taxonomy, namespaces and evidence codes in each of the given swissprot .dat files

    Returns a dictionary mapping filenames (keys) to protein counts (values)
    for the given taxonomy.
    """
    taxons = [str(taxon_id) for taxon_id in taxons]
    sp_count_by_namespace = []

    for filename in sp_file_list:
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
    #namespaces = list(data.get(list(data)[0]))
    namespaces = set([k[1] for k in [list(s.get("counts").keys()) for s in data][0]])
    print(namespaces)

    preheader = f"| COUNTS OF PROTEINS WITH ANNOTATION {' ' * (col1_width+10)}|"
    table_header = "| FILE{} | TAXON | {}    | {}    | {}    | GROWTH |"
    table_row_sep = (
        f"{'-' * col1_width}------------------------------------------------"
    )
    table_row_pattern = "| {} | {}  | {:6,} | {:6,} | {:6,} | {} |"

    print(table_row_sep)
    print(preheader)
    print(table_row_sep)
    print(table_header.format(" " * (col1_width - 4), *namespaces))
    print(table_row_sep)
    previous_row_sum = None

    counts = [row.get("counts") for row in data]

    #for key, counts in data.items():

    print(counts)

    #for key, row_count in counts.items():
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
