"""
Enzyme-to-Reaction Mapper
Extracts enzyme annotations from SBML models and maps them to reaction IDs.
"""

import libsbml
import pandas as pd
from typing import Dict, List, Set
import re

def extract_enzyme_annotations(xml_file_path: str) -> Dict[str, List[str]]:
    """
    Extract enzyme annotations from SBML model.
    
    Returns:
        Dictionary mapping enzyme names/EC numbers to list of reaction IDs
        Format: {enzyme_name: [reaction_id1, reaction_id2, ...]}
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(xml_file_path)
    model = document.getModel()
    
    if model is None:
        raise ValueError(f"Could not parse {xml_file_path}")
    
    enzyme_to_reactions = {}
    
    # Extract annotations from reactions
    for i in range(model.getNumReactions()):
        reaction = model.getReaction(i)
        reaction_id = reaction.getId()
        
        # Check for enzyme annotations in notes/annotations
        notes = reaction.getNotesString()
        annotation = reaction.getAnnotationString()
        
        # Look for EC numbers (EC:1.2.3.4 format)
        ec_numbers = re.findall(r'EC:\s*(\d+\.\d+\.\d+\.\d+)', notes + annotation)
        
        # Look for enzyme names in notes
        # Common patterns: "catalyzed by X", "enzyme: X", etc.
        enzyme_patterns = [
            r'enzyme[:\s]+([A-Za-z0-9\s\-]+)',
            r'catalyzed by[:\s]+([A-Za-z0-9\s\-]+)',
            r'EC number[:\s]+(\d+\.\d+\.\d+\.\d+)'
        ]
        
        enzymes_found = set()
        
        # Extract EC numbers
        for ec in ec_numbers:
            enzymes_found.add(f"EC:{ec}")
            if f"EC:{ec}" not in enzyme_to_reactions:
                enzyme_to_reactions[f"EC:{ec}"] = []
            enzyme_to_reactions[f"EC:{ec}"].append(reaction_id)
        
        # Extract enzyme names from patterns
        for pattern in enzyme_patterns:
            matches = re.findall(pattern, notes + annotation, re.IGNORECASE)
            for match in matches:
                enzyme_name = match.strip()
                if enzyme_name and len(enzyme_name) > 2:
                    enzymes_found.add(enzyme_name)
                    if enzyme_name not in enzyme_to_reactions:
                        enzyme_to_reactions[enzyme_name] = []
                    enzyme_to_reactions[enzyme_name].append(reaction_id)
        
        # Also check reaction name - sometimes enzyme name is in reaction name
        reaction_name = reaction.getName()
        if reaction_name:
            # Look for common enzyme name patterns in reaction names
            enzyme_in_name = re.findall(r'([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s+(?:synthase|reductase|dehydrogenase|kinase|phosphatase)', reaction_name)
            for enzyme in enzyme_in_name:
                if enzyme not in enzyme_to_reactions:
                    enzyme_to_reactions[enzyme] = []
                enzyme_to_reactions[enzyme].append(reaction_id)
    
    return enzyme_to_reactions


def create_enzyme_reaction_map(model_files: List[str]) -> pd.DataFrame:
    """
    Create a comprehensive mapping of enzymes to reactions across all models.
    
    Args:
        model_files: List of SBML file paths
        
    Returns:
        DataFrame with columns: [species, enzyme_name, reaction_id, reaction_name]
    """
    all_mappings = []
    
    species_map = {
        'M. gottschalkii.xml': 'MGK',
        'P. ruminicola.xml': 'PRM',
        'R. flavefaciens.xml': 'RFL'
    }
    
    for model_file in model_files:
        species = species_map.get(model_file, 'UNKNOWN')
        enzyme_reactions = extract_enzyme_annotations(f'model files/{model_file}')
        
        # Also get reaction names
        reader = libsbml.SBMLReader()
        document = reader.readSBML(f'model files/{model_file}')
        model = document.getModel()
        reaction_names = {r.getId(): r.getName() for r in [model.getReaction(i) for i in range(model.getNumReactions())]}
        
        for enzyme, reaction_ids in enzyme_reactions.items():
            for reaction_id in reaction_ids:
                reaction_name = reaction_names.get(reaction_id, '')
                all_mappings.append({
                    'species': species,
                    'enzyme_name': enzyme,
                    'reaction_id': reaction_id,
                    'reaction_name': reaction_name
                })
    
    return pd.DataFrame(all_mappings)


def map_docking_to_reactions(docking_df: pd.DataFrame, 
                             enzyme_reaction_map: pd.DataFrame) -> pd.DataFrame:
    """
    Map docking results (which have enzyme names) to reaction IDs.
    
    Args:
        docking_df: DataFrame with docking results (must have 'enzyme_name' column)
        enzyme_reaction_map: DataFrame from create_enzyme_reaction_map()
        
    Returns:
        DataFrame with docking results mapped to reaction IDs
    """
    # Merge on enzyme name
    merged = docking_df.merge(
        enzyme_reaction_map,
        on=['enzyme_name', 'species'],
        how='left'
    )
    
    # Handle cases where enzyme name doesn't match exactly
    # Try fuzzy matching or manual mapping if needed
    
    return merged


if __name__ == "__main__":
    # Example usage
    model_files = ['M. gottschalkii.xml', 'P. ruminicola.xml', 'R. flavefaciens.xml']
    
    print("Creating enzyme-to-reaction mapping...")
    enzyme_map = create_enzyme_reaction_map(model_files)
    
    # Save to CSV
    enzyme_map.to_csv('docking_data/enzyme_reaction_map.csv', index=False)
    print(f"Saved mapping with {len(enzyme_map)} enzyme-reaction pairs")
    print("\nSample mappings:")
    print(enzyme_map.head(20))

