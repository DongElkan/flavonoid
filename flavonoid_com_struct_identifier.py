"""
Complex structure identifier

This module identifies complex flavonoid skeleton, including
    *   biflavonoids
    *   theaflavins
    *   chalcone polymers
    *   spirobiflavonoids
"""


def _chalcone_polymers(molobj, cchains):
    """
    Identify chalcone polymers generated by Diels-Alder Reaction
    """
    # Use substructure search to search complex structure for
    # identifying chalcone polymers.
    pass