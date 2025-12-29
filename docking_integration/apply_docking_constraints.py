"""
Apply Docking-Based Constraints to the Cowmunity Model
This module integrates with the existing Cowmunity.py Variables() function.
"""

from .parse_docking_results import create_constraint_dict
from typing import Dict, Optional

def apply_docking_treatment(v_mgk, v_prm, v_rfl, 
                           molecule_name: str,
                           docking_csv_path: str,
                           enzyme_map_path: Optional[str] = None,
                           methane: str = 'variable'):
    """
    Apply docking-based constraints to model variables.
    
    Args:
        v_mgk, v_prm, v_rfl: GAMSpy Variable objects for each species
        molecule_name: Name of the drug molecule
        docking_csv_path: Path to CSV file with docking results
        enzyme_map_path: Optional path to enzyme-reaction mapping CSV
        methane: 'variable' or 'fixed' methane mode
    
    This function modifies the variable bounds based on docking data.
    """
    print(f"Applying docking-based constraints for {molecule_name}...")
    
    # Get constraints from docking data
    constraints = create_constraint_dict(molecule_name, docking_csv_path, enzyme_map_path)
    
    # Apply constraints to each species
    species_vars = {
        'MGK': v_mgk,
        'PRM': v_prm,
        'RFL': v_rfl
    }
    
    for species, reaction_constraints in constraints.items():
        var = species_vars[species]
        print(f"\nApplying constraints to {species}:")
        
        for reaction_id, constraint_data in reaction_constraints.items():
            inhibition_factor = constraint_data['inhibition_factor']
            inhibition_type = constraint_data['inhibition_type']
            
            # Check if reaction exists in the model
            try:
                current_lb = var.lo[reaction_id] if hasattr(var, 'lo') else None
                current_ub = var.up[reaction_id] if hasattr(var, 'up') else None
                
                if current_lb is None or current_ub is None:
                    print(f"  Warning: Reaction {reaction_id} not found in {species} model")
                    continue
                
                # Apply inhibition based on type
                if inhibition_type.lower() in ['competitive', 'competitive inhibitor']:
                    # Competitive: reduce upper bound
                    try:
                        current_ub_val = float(current_ub) if hasattr(current_ub, '__float__') else current_ub
                        new_ub = current_ub_val * inhibition_factor
                        var.up[reaction_id] = new_ub
                        print(f"  {reaction_id}: UB {current_ub_val:.4f} → {new_ub:.4f} (competitive)")
                    except Exception as e:
                        # If we can't get numeric value, just set it
                        var.up[reaction_id] = var.up[reaction_id] * inhibition_factor
                        print(f"  {reaction_id}: UB reduced by {(1-inhibition_factor)*100:.1f}% (competitive)")
                
                elif inhibition_type.lower() in ['non-competitive', 'noncompetitive', 'non competitive']:
                    # Non-competitive: reduce both bounds
                    try:
                        current_lb_val = float(current_lb) if hasattr(current_lb, '__float__') else current_lb
                        current_ub_val = float(current_ub) if hasattr(current_ub, '__float__') else current_ub
                        new_lb = current_lb_val * inhibition_factor
                        new_ub = current_ub_val * inhibition_factor
                        var.lo[reaction_id] = new_lb
                        var.up[reaction_id] = new_ub
                        print(f"  {reaction_id}: LB {current_lb_val:.4f} → {new_lb:.4f}, UB {current_ub_val:.4f} → {new_ub:.4f} (non-competitive)")
                    except Exception as e:
                        var.lo[reaction_id] = var.lo[reaction_id] * inhibition_factor
                        var.up[reaction_id] = var.up[reaction_id] * inhibition_factor
                        print(f"  {reaction_id}: Bounds reduced by {(1-inhibition_factor)*100:.1f}% (non-competitive)")
                
                elif inhibition_type.lower() in ['uncompetitive', 'un competitive']:
                    # Uncompetitive: reduce both bounds proportionally
                    new_lb = current_lb * inhibition_factor
                    new_ub = current_ub * inhibition_factor
                    var.lo[reaction_id] = new_lb
                    var.up[reaction_id] = new_ub
                    print(f"  {reaction_id}: LB {current_lb:.4f} → {new_lb:.4f}, UB {current_ub:.4f} → {new_ub:.4f} (uncompetitive)")
                
                elif inhibition_type.lower() in ['activator', 'activation']:
                    # Activator: increase bounds
                    activation_factor = inhibition_factor  # > 1.0
                    new_lb = current_lb * activation_factor if current_lb < 0 else current_lb
                    new_ub = current_ub * activation_factor
                    var.lo[reaction_id] = new_lb
                    var.up[reaction_id] = new_ub
                    print(f"  {reaction_id}: UB {current_ub:.4f} → {new_ub:.4f} (activator)")
                
                else:
                    # Default: treat as competitive
                    new_ub = current_ub * inhibition_factor
                    var.up[reaction_id] = new_ub
                    print(f"  {reaction_id}: UB {current_ub:.4f} → {new_ub:.4f} (default)")
            
            except KeyError:
                print(f"  Warning: Reaction {reaction_id} not found in {species} model")
                continue
            except Exception as e:
                print(f"  Error applying constraint to {reaction_id}: {e}")
                continue
    
    print(f"\nDocking-based constraints applied for {molecule_name}")


# Integration function to be called from Cowmunity.py
def apply_docking_treatment_integrated(treatment: str, 
                                       v_mgk, v_prm, v_rfl,
                                       methane: str = 'variable'):
    """
    Integration function that maps treatment names to docking CSV files.
    Modify this to point to your actual docking data files.
    """
    docking_files = {
        'imidazole': 'docking_data/imidazole_docking.csv',
        'l-carnitine': 'docking_data/l-carnitine_docking.csv',
        'methyl jasmonate': 'docking_data/methyl_jasmonate_docking.csv',
        'propylpyrazine': 'docking_data/propylpyrazine_docking.csv'
    }
    
    if treatment in docking_files:
        docking_path = docking_files[treatment]
        apply_docking_treatment(
            v_mgk, v_prm, v_rfl,
            molecule_name=treatment,
            docking_csv_path=docking_path,
            methane=methane
        )
    else:
        print(f"No docking data available for treatment: {treatment}")

