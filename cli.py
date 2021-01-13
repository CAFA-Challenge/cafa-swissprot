""" Functions for the command-line entrypoints defined in setup.py
for the CAFA sprot targets module.

Three entrypoints are defined:
CAFA_experimental_growth
CAFA_print_annotation_counts
CAFA_generate_no_exp_files

"""
import os
from pathlib import Path
import yaml
import click
from experimental_growth import (
    count_annotated_proteins_from_files,
    print_annotation_counts_table,
)
from filter_sprot_species import (
    species_filter,
    filter_sprot_by_taxonomies,
    generate_protein_ids_mapping,
)
from utils import TAXONOMY_LOOKUP


@click.command()
@click.argument("config_handle", type=click.File("r"))
def experimental_growth(config_handle):
    """Collects counts of annotations for the given taxonomy per ontology namespace
    filtered by the provided evidence codes. The counts are printed to stdout.
    """
    conf = yaml.load(config_handle, Loader=yaml.FullLoader)
    namespaces = conf.get("ontologies")
    evidence_codes = conf.get("evidence_codes")
    swissprot_file_list = conf.get("input_files")
    taxons = conf.get("taxonomies")
    if taxons is None:
        taxon = conf.get("taxonomy")

        if taxon is None:
            raise Exception("either taxon or taxonomies is required")

        taxons = (taxon,)

    growth_counts = count_annotated_proteins_from_files(
        swissprot_file_list,
        taxons,
        allowed_namespaces=namespaces,
        allowed_evidence_codes=evidence_codes,
    )

    keys = ("filename", "taxon_id", "taxon_name") + tuple(namespaces)
    print("\t".join(keys))

    # Process the results:
    # 1. Derive nested lists from the nested dicts
    # 2. Derive SORTED nested lists from the data in step #1
    _unsorted = []

    for count_record in growth_counts:
        sprot_filename = count_record.get("filename")
        for key, counts in count_record.get("counts").items():
            taxon, namespace = key
            _unsorted.append([sprot_filename, taxon, namespace, counts])

    # Loop through the configuration taxons, filelist and ontologies in order
    # to main the original order of things:
    for taxon in taxons:
        for sprot_filename in swissprot_file_list:

            counts = [
                row
                for row in _unsorted
                if row[0] == sprot_filename and row[1] == str(taxon)
            ]

            # create a new list merging the individual namespace counts into a
            # single row
            combined_row = [sprot_filename, str(taxon), TAXONOMY_LOOKUP.get(str(taxon))]

            for ontology in namespaces:
                # row[-1] is the protein count, row[-2] is the ontology namespace
                combined_row += [str(row[-1]) for row in counts if row[-2] == ontology]

            print("\t".join(combined_row))


@click.command()
@click.argument("config_handle", type=click.File("r"))
def print_annotation_counts(config_handle):
    """Collects counts of annotations for the given taxonomy per ontology namespace
    filtered by the provided evidence codes. The counts are printed to stdout.
    """
    conf = yaml.load(config_handle, Loader=yaml.FullLoader)
    taxons = conf.get("taxonomies") or (conf.get("taxon"),)
    namespaces = conf.get("ontologies")
    evidence_codes = conf.get("evidence_codes")
    swissprot_file_list = conf.get("input_files")

    data = count_annotated_proteins_from_files(
        swissprot_file_list,
        taxons,
        allowed_namespaces=namespaces,
        allowed_evidence_codes=evidence_codes,
    )
    print_annotation_counts_table(data, taxons)


@click.command()
@click.argument("config_handle", type=click.File("r"))
@click.option(
    "-q",
    "--quiet",
    "quiet",
    is_flag=True,
    default=False,
    help="suppress stdout progress messages",
)
def generate_no_exp_files(config_handle, quiet=False):
    """Parses a swissprot file and extracts protein data for the given taxonomy
    where there is little/no experimental annotation for the relevant GO namespaces
    and evidence codes. That protein data is written to a series of fasta files.

    CONFIG_HANDLE is the filepath for a yaml file containing various parameters. See
    filter_sprot_species_example.yml in the git repo for more info.
    """

    def _print(message: str):
        """ simple wrapper around print() obeying the 'quiet' flag """
        if not quiet:
            print(message)

    conf = yaml.load(config_handle, Loader=yaml.FullLoader)

    taxonomies = conf.get("taxonomies")

    if taxonomies is None:
        taxon = conf.get("taxonomy")
        if taxon is None:
            raise Exception("either taxon or taxonomies is required")
        else:
            taxonomies = (taxon,)

    namespaces = conf.get("ontologies")
    allowed_evidence_codes = conf.get("allowed_evidence_codes")
    sprot_file = conf.get("sprot_file")
    output_directory = conf.get("output_directory")

    _print(f"Parsing {sprot_file} for taxonomies {taxonomies}")

    output_dir_as_path = Path(output_directory)
    if not output_dir_as_path.is_dir():
        _print(f"\tMaking output directory {output_directory}")
        os.makedirs(output_dir_as_path, exist_ok=True)

    _print(f"Opening {sprot_file}")
    with open(sprot_file, "r") as sprot_handle:
        #species_filter(sprot_handle, taxonomies=taxonomies)
        _print("Filtering by GO Namespace and Evidence Code")
        if namespaces:
            _print(f"\tUsing GO namespaces: {namespaces}")
        if allowed_evidence_codes:
            _print(f"\tUsing ALLOWED evidence codes: {allowed_evidence_codes}")

        _print(f"Writing results to {output_directory}")
        filter_sprot_by_taxonomies(
            sprot_handle,
            output_dir=output_directory,
            taxonomies=taxonomies,
            namespaces=namespaces,
            allowed_evidence_codes=allowed_evidence_codes,
        )
        _print(f"Generating Swissprot ID => CAFA ID map file(s) for {output_directory}")
        generate_protein_ids_mapping(taxonomies, output_directory)
        _print("")
