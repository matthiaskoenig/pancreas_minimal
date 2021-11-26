"""Model annotations."""
from sbmlutils.metadata import *

compartments = {
    # pancreas
    "pa": [
        (BQB.IS, "fma/FMA:7198"),
        (BQB.IS, "ncit/C12393"),
        (BQB.IS, "uberon/UBERON:0001264"),
        (BQB.IS, "bto/BTO:0000988"),
    ],
    "plasma": [
        (BQB.IS, "ncit/C13356"),
        (BQB.IS, "bto/BTO:0000131"),
    ],
    "blood": [
        (BQB.IS, "fma/FMA:62970"),
        (BQB.IS, "ncit/C12434"),
        (BQB.IS, "bto/BTO:0000089"),
    ],
    "plasma membrane": [
        (BQB.IS, "fma/FMA:63841"),
        (BQB.IS, "ncit/C13735"),
    ],
}

species = {
    "glc": [
        (BQB.IS, "chebi/CHEBI:17234"),
        (BQB.IS, "ncit/C2831"),
    ],
    "lac": [
        (BQB.IS, "chebi/CHEBI:24996"),
    ],
    "ins": [
        (BQB.IS, "chebi/CHEBI:24996"),
        (BQB.IS, "ncit/C2271"),
        (BQB.IS, "uniprot/P01308"),
    ],
    "cpep": [
        (BQB.IS, "chebi/CHEBI:80332"),
        (BQB.IS, "ncit/C94608"),
        (BQB.IS_VERSION_OF, "uniprot/P01308"),
    ],
}
