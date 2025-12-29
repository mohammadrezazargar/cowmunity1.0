"""
Docking Integration Module for Cowmunity Model

This module integrates molecular docking data into metabolic flux constraints.
"""

from .parse_docking_results import (
    parse_docking_csv,
    calculate_inhibition_factor,
    create_constraint_dict,
    group_by_reaction
)

from .apply_docking_constraints import (
    apply_docking_treatment,
    apply_docking_treatment_integrated
)

from .enzyme_reaction_mapper import (
    extract_enzyme_annotations,
    create_enzyme_reaction_map,
    map_docking_to_reactions
)

__all__ = [
    'parse_docking_csv',
    'calculate_inhibition_factor',
    'create_constraint_dict',
    'group_by_reaction',
    'apply_docking_treatment',
    'apply_docking_treatment_integrated',
    'extract_enzyme_annotations',
    'create_enzyme_reaction_map',
    'map_docking_to_reactions'
]

