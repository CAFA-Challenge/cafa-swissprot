import sys
import os
from pathlib import Path
from typing import Iterable, Optional, TextIO
from operator import itemgetter
from itertools import count
from Bio import SeqIO, SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from utils import (
    has_go_annotation,
    parse_go,
    NAMESPACES,
    ALLOWED_EVIDENCE_CODES,
)


def go_evidence_code_filter(
    sprot_record,
    namespaces: Optional[Iterable] = NAMESPACES,
    # namespace: str,
    allowed_evidence_codes: Optional[Iterable] = ALLOWED_EVIDENCE_CODES,
    use_has_go: bool = True,
) -> dict:
    """
    Returns a dict mapping ontology namespaces to booleans
    1. False values when there is no go association for an allowed namespace
    2. True values otherwise
    """
    in_allowed_ontologies = {namespace: True for namespace in namespaces}

    if use_has_go and not has_go_annotation(sprot_record):
        return in_allowed_ontologies

    go_list = parse_go(sprot_record)

    # Check that all GO experimental annotations are in allowed namespaces:
    for namespace in in_allowed_ontologies.keys():
        namespace_go_list = [
            term.evidence_code in allowed_evidence_codes
            for term in go_list
            if term.ontology == namespace
        ]
        in_allowed_ontologies[namespace] = all(namespace_go_list)

    return in_allowed_ontologies


def species_filter(
    sp_handle: TextIO,
    #output_handle: TextIO = None,
    output_directory: Optional[str] = ".",
    taxonomies=Iterable[str],
):
    """Filters swissprot data on the given taxonomy ID and
    writes the filtered data either to the specified file-like object
    """
    taxonomies = {str(taxon) for taxon in taxonomies}
    ids = get_id_sequence()

    output_dir_as_path = Path(output_directory)
    species_filter_output_path = "sp_species.{}.tfa"

    write_handles = {
        taxon: open(output_dir_as_path / species_filter_output_path.format(taxon), "w")
        for taxon in taxonomies
    }

    '''
    # Write results either to a (1) specified file-like handle or to (2) stdout:
    if output_handle is not None:
        write_handle = output_handle
    else:
        write_handle = sys.stdout
    '''

    sp_handle.seek(0)

    # filter the swissprot data by taxonomy:
    sprot_records = [
        record
        for record in SwissProt.parse(sp_handle)
        if len(taxonomies.intersection(set(record.taxonomy_id))) > 0
    ]

    for record in sprot_records:
        # if taxon not in record.taxonomy_id:
        #   continue
        tax_intersect = taxonomies.intersection(set(record.taxonomy_id))
        # Is this valid? Can we assume there's only one taxon here?
        _taxon = list(tax_intersect)[0]

        protein_seq = SeqRecord(
            Seq(record.sequence),
            id=f"T{_taxon}{next(ids):0>7}",
            description=record.entry_name,
        )
        SeqIO.write([protein_seq], write_handles[_taxon], "fasta")


def get_id_sequence():
    """ generator for seq IDs """
    for i in count(start=1):
        yield i


def write_no_exp_files(
    sprot_handle: TextIO,
    taxonomy: str,
    namespaces: Iterable[str],
    output_dir: Optional[str],
    allowed_evidence_codes: Optional[str],
):

    """filters Swiss Prot data for a single TAXON as well as ontology
    namespaces and allowed evidence codes.

    Writes 1 or more (based on number of provided namespaces) output files
    """
    output_dir_path = Path(output_dir or ".")

    taxonomy = str(taxonomy)
    # instance of the ID generator function:
    ids = get_id_sequence()
    write_filename_pattern = "sp_species.{}.{}.noexp.tfa"
    out_handles = {
        "ALL": open(
            output_dir_path / Path(write_filename_pattern.format(taxonomy, "all")), "w"
        )
    }

    for namespace in namespaces:
        out_handles[namespace] = open(
            output_dir_path / Path(write_filename_pattern.format(taxonomy, namespace)),
            "w",
        )

    sprot_handle.seek(0)
    # Filter Swiss Prot data by taxon:
    proteins = [
        record
        for record in SwissProt.parse(sprot_handle)
        if taxonomy in record.taxonomy_id
    ]

    for record in proteins:
        assert taxonomy in record.taxonomy_id
        # This gives us the next integer in the sequence:
        _id = next(ids)

        # go_evidence_code_filter() returns a dict mapping ontologies (keys) to boolean values
        is_in_allowed = go_evidence_code_filter(
            record, namespaces=namespaces, allowed_evidence_codes=allowed_evidence_codes
        )

        for ontology, is_allowed in is_in_allowed.items():
            if is_allowed is False:
                continue

            protein_seq = SeqRecord(
                Seq(record.sequence),
                id=f"T{taxonomy}{_id:0>7}",
                description=record.entry_name,
            )
            SeqIO.write(protein_seq, out_handles[ontology], "fasta")

        if all(is_in_allowed.values()):
            SeqIO.write(protein_seq, out_handles["ALL"], "fasta")

    for handle in out_handles.values():
        handle.close()


def filter_sprot_by_taxonomies(
    sp_handle: TextIO,
    taxonomies: Iterable[str],
    output_dir: Optional[str] = ".",
    namespaces: Optional[Iterable] = NAMESPACES,
    allowed_evidence_codes: Optional[Iterable] = ALLOWED_EVIDENCE_CODES,
):
    """ Writes 4 output files bases on GO namespaces """

    for taxonomy in taxonomies:
        write_no_exp_files(
            sprot_handle=sp_handle,
            taxonomy=str(taxonomy),
            namespaces=namespaces,
            output_dir=output_dir,
            allowed_evidence_codes=allowed_evidence_codes,
        )


# def generate_protein_ids_mapping(taxon_id: str, read_directory: str):
def generate_protein_ids_mapping(taxonomies: Iterable[str], read_directory: str):
    """Writes a .txt file that contains a mapping of generated CAFA protein IDs
    to the original protein IDs from the swissprot .dat file(s).
    """
    file_extensions_of_interest = (".tfa", ".fa", ".fasta")
    taxonomies = {str(taxon) for taxon in taxonomies}

    for root, _, filenames in os.walk(read_directory):
        for taxon in taxonomies:
            collected_id_tuples = []
            taxon_files = [
                Path(root, filename)
                for filename in filenames
                if taxon in filename
                and Path(root, filename).suffix in file_extensions_of_interest
            ]

            for filepath in taxon_files:

                with open(filepath, "r") as read_handle:
                    lines = [line for line in read_handle if line.startswith(">")]
                    # lines of interest have this form:
                    # >T96060000024 A20A1_HUMAN
                    for line in lines:
                        cafa_id, sprot_id = line.rstrip()[1:].split()
                        collected_id_tuples.append((cafa_id, sprot_id))

            write_dir_path = Path(read_directory)
            with open(
                write_dir_path / Path(f"{taxon}_id_map.txt"), "w"
            ) as write_handle:
                # Make sure the ID pairs are unique and sorted by the generated CAFA ID:
                collected_id_tuples = sorted(
                    tuple(set(collected_id_tuples)), key=itemgetter(0)
                )

                for cafa_id, sprot_id in collected_id_tuples:
                    write_handle.write(f"{cafa_id}\t{sprot_id}\n")
