from typing import Optional, Iterable
from collections import namedtuple


# the data type for storing protein go annotations:
Annotation = namedtuple("Annotation", ("protein", "go_id", "ontology", "evidence_code"))

namespace_lookup = {"P": "BPO", "C": "CCO", "F": "MFO"}
NAMESPACES = [v for k, v in namespace_lookup.items()]

EXPERIMENTAL_EVIDENCE_CODES = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")

ALLOWED_EVIDENCE_CODES = (
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

# Taxon codes can be found at https://www.uniprot.org/docs/speclist.txt
TAXONOMY_LOOKUP = {
    "9606": "Homo sapiens",
    "10090": "Mus musculus",
    "10116": "Rattus norvegicus",
    "3702": "Arabidopsis thaliana",
    "83333": "Escherichia coli",
    "7227": "Drosophila melanogaster",
    "287": "Pseudomonas aeruginosa",
    "559292": "Saccharomyces cerevisiae",
    "284812": "Schizosaccharomyces pombe",
    "7955": "Danio rerio",
    "44689": "Dictyostelium discoideum",
    "243273": "Mycoplasma genitalium",
    "6239": "Caenorhabditis elegans",
    "226900": "Bacillus cereus",
    "4577": "Zea Mays",
    "9823": "Sus scrofa",
    "99287": "Salmonella typhymurium",
}


def has_go_annotation(sp_record) -> bool:
    """ returns True/False based on existence of GO reference in the swissprot data """
    for xref in sp_record.cross_references:
        if xref[0] == "GO":
            return True

    return False


def parse_go(
    sp_rec,
    namespace_filter: Optional[Iterable] = None,
    evidence_code_filter: Optional[Iterable] = None,
) -> list:
    """Given a swissprot record object, returns a list of GO terms (as dicts)
    taken from the swissprot record object.cross_references data.
    """
    if namespace_filter is not None:
        namespace_filter = [namespace.upper() for namespace in namespace_filter]

    annotations = []
    go_references = (ref for ref in sp_rec.cross_references if ref[0] == "GO")

    for ref in go_references:
        # a ref is a tuple with this form:
        # ('GO', 'GO:0015629', 'C:actin cytoskeleton', 'IDA:BHF-UCL'),
        go_id, ontology_namespace, evidence_code = ref[1:4]
        ontology_namespace = namespace_lookup.get(ontology_namespace.split(":")[0])
        evidence_code = evidence_code.split(":")[0]

        if namespace_filter is not None and ontology_namespace not in namespace_filter:
            continue

        if (
            evidence_code_filter is not None
            and evidence_code not in evidence_code_filter
        ):
            continue

        annotations.append(
            Annotation(sp_rec.entry_name, go_id, ontology_namespace, evidence_code)
            # {"go_id": go_id, "ontology": ontology_namespace, "evidence": evidence_code}
        )

    return annotations
