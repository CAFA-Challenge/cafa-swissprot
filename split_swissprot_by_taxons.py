from typing import Iterable
from pathlib import Path
from Bio import SeqIO


def split_swissprot_by_taxons(
    swissprot_filepath: str, taxons: Iterable[str], write_dir_path="./"
) -> None:
    """Reads a given Swissprot file and writes N (one per taxon)
    Swissprot-formatted files based on the given taxons IDs.
    """
    taxons = {str(species_id) for species_id in taxons}
    write_path = Path(write_dir_path)

    write_handles = {
        taxon_id: open(write_path.joinpath(f"swissprot_{taxon_id}.dat"), "w")
        for taxon_id in taxons
    }

    # This is expensive b/c it loads the entire file into memory,
    # but it is the easiest way to output the data in swissprot format
    swissprot_dict = SeqIO.index(swissprot_filepath, "swiss")

    for key, val in swissprot_dict.items():
        taxon_id = val.__dict__.get("annotations").get("ncbi_taxid")
        _intersect = set(taxon_id).intersection(taxons)

        if len(_intersect) > 0:
            _taxon = list(_intersect)[0]
            write_handles[_taxon].write(swissprot_dict.get_raw(key).decode())

    for _, file_handle in write_handles.items():
        file_handle.close()


if __name__ == "__main__":
    SPROT_FILEPATH = (
        "/home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot.dat"
    )
    TAXONS_OF_INTEREST = (
        9606,  # human
        4577,  # maize
    )

    split_swissprot_by_taxons(
        swissprot_filepath=SPROT_FILEPATH,
        taxons=TAXONS_OF_INTEREST,
        write_dir_path="TRASH",
    )
