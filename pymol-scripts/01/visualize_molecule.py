#!/usr/bin/env python
import sys
import os
import argparse
import pymol
from pymol import cmd, stored

def find_free_cysteines(obj_name, distance_cutoff=2.5):
    """
    Find free cysteines (not in disulfide bonds)
    
    Args:
        obj_name: PyMOL object name
        distance_cutoff: Maximum S-S distance for disulfide bond (default 2.5 Å)
    
    Returns:
        List of (chain, resi, min_distance) tuples
    """
    free_cys = []
    cmd.select("all_cys", f"{obj_name} and resn CYS")
    stored.cys_residues = []
    cmd.iterate("all_cys and name CA", "stored.cys_residues.append((chain, resi))")
    
    for chain, resi in stored.cys_residues:
        sel = f"{obj_name} and chain {chain} and resi {resi} and name SG"
        
        # Get SG coordinates
        stored.sg_coords = []
        cmd.iterate_state(1, sel, "stored.sg_coords.append([x, y, z])")
        
        if not stored.sg_coords:
            continue
            
        # Find nearest SG atom
        min_dist = float('inf')
        stored.other_sg = []
        cmd.iterate_state(1, f"{obj_name} and resn CYS and name SG and not ({sel})", 
                         "stored.other_sg.append([x, y, z])")
        
        for other_coords in stored.other_sg:
            dist = sum((stored.sg_coords[0][i] - other_coords[i])**2 for i in range(3))**0.5
            min_dist = min(min_dist, dist)
        
        # If no nearby SG or distance > cutoff, it's free
        if min_dist > distance_cutoff:
            free_cys.append((chain, resi, min_dist if min_dist != float('inf') else None))
    
    cmd.delete("all_cys")
    return free_cys

def visualize_molecule(pdb_file, output_dir=None, show_sidechain=False, quality='low', view_style='cartoon'):
    if not os.path.exists(pdb_file):
        print(f"Error: File {pdb_file} not found")
        return
    
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    if output_dir is None:
        output_dir = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Set quality parameters
    quality_settings = {
        'low': {'width': 600, 'height': 450, 'dpi': 100, 'ray': False},
        'medium': {'width': 800, 'height': 600, 'dpi': 150, 'ray': True},
        'high': {'width': 1200, 'height': 900, 'dpi': 200, 'ray': True}
    }
    settings = quality_settings.get(quality, quality_settings['low'])
    
    cmd.reinitialize()
    cmd.load(pdb_file, "molecule")
    cmd.remove("hydro")
    
    # Color chains with Mirabo-inspired palette (high contrast)
    colors = [
        [0.00, 0.23, 0.64],  # #003BA3 - Deep blue (Mirabo primary)
        [0.20, 0.80, 0.20],  # #32CD32 - Lime green (Mirabo accent)
        [1.00, 0.50, 0.00],  # Orange
        [0.60, 0.20, 0.80],  # Purple
        [0.00, 0.80, 0.80],  # Cyan
        [1.00, 0.75, 0.00],  # Gold
        [0.80, 0.00, 0.40],  # Magenta
        [0.29, 0.55, 1.00],  # #4B8DFF - Light blue (Mirabo secondary)
    ]
    stored.chains = []
    cmd.iterate("molecule and name CA", "stored.chains.append(chain)")
    unique_chains = list(dict.fromkeys(stored.chains))  # Preserve order
    
    for i, chain in enumerate(unique_chains):
        color_name = f"mirabo_color_{i}"
        cmd.set_color(color_name, colors[i % len(colors)])
        cmd.color(color_name, f"molecule and chain {chain}")
    
    # Apply view style
    cmd.show(view_style, "molecule")
    
    # Mark cysteines
    cmd.select("all_cysteines", "molecule and resn CYS")
    cmd.show("sticks", "all_cysteines")
    cmd.color("yellow", "all_cysteines")
    
    # Find and mark free cysteines
    free_cys = find_free_cysteines("molecule", distance_cutoff=2.5)
    if free_cys:
        for chain, resi, min_dist in free_cys:
            cmd.select("free_cys", f"molecule and chain {chain} and resi {resi}")
            cmd.show("sticks", "free_cys")
            cmd.color("red", "free_cys")
            cmd.set("stick_radius", 0.3, "free_cys")
    
    if show_sidechain:
        cmd.show("lines", "molecule and sidechain")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    
    # Standard views
    cmd.orient("molecule")
    cmd.png(f"{output_dir}/{base_name}_front.png", width=1920, height=1080, dpi=300, ray=1)
    
    cmd.turn("x", 90)
    cmd.png(f"{output_dir}/{base_name}_top.png", width=1920, height=1080, dpi=300, ray=1)
    
    cmd.orient("molecule")
    cmd.turn("y", 90)
    cmd.png(f"{output_dir}/{base_name}_side.png", width=1920, height=1080, dpi=300, ray=1)
    
    cmd.orient("molecule")
    cmd.zoom("molecule", 10)
    cmd.png(f"{output_dir}/{base_name}_overview.png", width=1920, height=1080, dpi=300, ray=1)
    
    # Free cysteine close-ups
    if free_cys:
        for idx, (chain, resi, min_dist) in enumerate(free_cys):
            cmd.orient(f"molecule and chain {chain} and resi {resi}")
            cmd.zoom(f"molecule and chain {chain} and resi {resi}", 15)
            cmd.label(f"molecule and chain {chain} and resi {resi} and name CA", f"'Chain {chain} CYS {resi}'")
            cmd.set("label_size", 20)
            cmd.set("label_color", "red")
            cmd.png(f"{output_dir}/{base_name}_free_cys_{idx+1}.png", width=1920, height=1080, dpi=300, ray=1)
            cmd.hide("labels")
    
    # Generate rotating GIF
    cmd.orient("molecule")
    cmd.zoom("molecule", 5)
    
    # Complete 360 degree rotation
    for i in range(36):
        cmd.turn("y", 10)
        cmd.rebuild()
        if settings['ray']:
            cmd.ray(settings['width'], settings['height'])
        cmd.png(f"{output_dir}/frame_{i:03d}.png", width=settings['width'], height=settings['height'], dpi=settings['dpi'])
    
    if free_cys:
        chain, resi, min_dist = free_cys[0]
        
        # Save the current view after rotation
        view_start = cmd.get_view()
        
        # Get the target view (zoomed to free cysteine)
        cmd.zoom(f"molecule and chain {chain} and resi {resi}", 2)
        view_end = cmd.get_view()
        
        # Interpolate between the two views
        total_frames = 40
        for i in range(total_frames):
            t = i / (total_frames - 1)  # 0 to 1
            # Linear interpolation for each view parameter
            view_current = tuple(
                view_start[j] + t * (view_end[j] - view_start[j])
                for j in range(len(view_start))
            )
            cmd.set_view(view_current)
            cmd.rebuild()
            if settings['ray']:
                cmd.ray(settings['width'], settings['height'])
            cmd.png(f"{output_dir}/frame_{36+i:03d}.png", width=settings['width'], height=settings['height'], dpi=settings['dpi'])
    
    gif_file = f"{output_dir}/{base_name}_rotation.gif"
    os.system(f"magick -dispose previous -delay 20 -loop 0 {output_dir}/frame_*.png {gif_file}")
    os.system(f"rm {output_dir}/frame_*.png")
    
    # Save PyMOL session for interactive viewing
    session_file = f"{output_dir}/{base_name}_session.pse"
    cmd.save(session_file)
    
    print(f"Visualization complete. Files saved to {output_dir}")
    print(f"Free cysteines found: {len(free_cys)}")
    if free_cys:
        for chain, resi, min_dist in free_cys:
            dist_str = f"{min_dist:.2f} Å" if min_dist else "N/A"
            print(f"  - Chain {chain}, Residue {resi} (nearest SG: {dist_str})")
    print(f"\nTo view interactively: pymol {session_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Visualize protein structures with PyMOL and identify free cysteines.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Description:
  This script loads a PDB file, visualizes the protein structure using PyMOL,
  and automatically identifies free (unpaired) cysteines that are not involved
  in disulfide bonds. It generates multiple views and an animated rotation.

View Styles:
  cartoon  - Secondary structure representation showing alpha helices and beta sheets
             Best for: Understanding overall protein fold and domain architecture
  
  surface  - Molecular surface representation showing the solvent-accessible surface
             Best for: Analyzing binding sites, protein-protein interfaces, cavities
  
  sticks   - Stick representation showing all atoms as connected cylinders
             Best for: Detailed examination of atomic positions and bonds
  
  spheres  - Space-filling spheres (CPK model) showing van der Waals radii
             Best for: Understanding molecular volume and packing
  
  lines    - Simple line representation connecting bonded atoms
             Best for: Quick overview with minimal visual clutter

Quality Settings:
  low      - 600x450 px, 100 DPI, no ray tracing (fast, smaller files)
  medium   - 800x600 px, 150 DPI, with ray tracing (balanced)
  high     - 1200x900 px, 200 DPI, with ray tracing (slow, publication quality)
  
  Note: Quality only affects the animated GIF. Static images are always 1920x1080 at 300 DPI.

Free Cysteine Detection:
  - Identifies cysteine residues not forming disulfide bonds (S-S distance > 2.5 Å)
  - Highlights free cysteines in RED with thicker sticks
  - Regular cysteines shown in YELLOW
  - Generates close-up views of each free cysteine
  - Reports chain, residue number, and nearest sulfur distance

Color Scheme:
  - Each protein chain is colored with a distinct Mirabo-inspired palette
  - High contrast colors for easy chain identification
  - Free cysteines: RED (potential conjugation sites)
  - All cysteines: YELLOW sticks overlaid on main structure

Output Files:
  <name>_front.png       - Front view of the molecule
  <name>_top.png         - Top view (90° rotation on X-axis)
  <name>_side.png        - Side view (90° rotation on Y-axis)
  <name>_overview.png    - Zoomed out overview
  <name>_free_cys_N.png  - Close-up of each free cysteine (if found)
  <name>_rotation.gif    - 360° rotation animation with zoom to first free cysteine
  <name>_session.pse     - PyMOL session file for interactive exploration

Examples:
  # Basic usage with default settings
  %(prog)s protein.pdb
  
  # Specify output directory
  %(prog)s protein.pdb output_images/
  
  # Surface view with high quality rendering
  %(prog)s protein.pdb --view surface --quality high
  
  # Show side chains with cartoon view
  %(prog)s protein.pdb --show-sidechain --view cartoon
  
  # Complete example with all options
  %(prog)s protein.pdb results/ --view surface --quality medium --show-sidechain

Tips:
  - Use 'cartoon' for presentations and overall structure analysis
  - Use 'surface' for binding site and interface analysis
  - Use 'sticks' or 'spheres' for detailed atomic-level examination
  - Open the .pse session file in PyMOL for interactive manipulation
  - Free cysteines are potential sites for antibody-drug conjugation (ADC)
        '''
    )
    
    parser.add_argument('pdb_file', 
                        metavar='PDB_FILE',
                        help='Path to input PDB file containing protein structure')
    
    parser.add_argument('output_dir', 
                        metavar='OUTPUT_DIR',
                        nargs='?', 
                        default=None,
                        help='Directory for output files (default: current directory)')
    
    parser.add_argument('--view', 
                        choices=['cartoon', 'surface', 'sticks', 'spheres', 'lines'],
                        default='cartoon',
                        metavar='STYLE',
                        help='Visualization style: cartoon (default), surface, sticks, spheres, or lines')
    
    parser.add_argument('--quality', 
                        choices=['low', 'medium', 'high'],
                        default='low',
                        metavar='LEVEL',
                        help='Rendering quality for animated GIF: low (default), medium, or high')
    
    parser.add_argument('--show-sidechain', 
                        action='store_true',
                        help='Display amino acid side chains as lines (adds detail to visualization)')
    
    args = parser.parse_args()
    
    pymol.finish_launching(['pymol', '-c'])
    visualize_molecule(args.pdb_file, args.output_dir, args.show_sidechain, args.quality, args.view)
