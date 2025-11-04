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

def mark_hydrophobic_hydrophilic(obj_name, mode='color'):
    """
    Mark hydrophobic and hydrophilic residues
    
    Args:
        obj_name: PyMOL object name
        mode: 'color', 'surface', or 'highlight'
    """
    hydrophobic = ['ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO']
    hydrophilic = ['SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU']
    
    # Define custom colors
    cmd.set_color('hydrophobic_color', [0.8, 0.6, 0.2])  # tan/ochre
    cmd.set_color('hydrophilic_color', [0.2, 0.5, 0.8])  # steel blue
    
    if mode == 'color':
        cmd.color('hydrophobic_color', f"{obj_name} and resn {'+'.join(hydrophobic)}")
        cmd.color('hydrophilic_color', f"{obj_name} and resn {'+'.join(hydrophilic)}")
    elif mode == 'surface':
        cmd.create('hydrophobic_surf', f"{obj_name} and resn {'+'.join(hydrophobic)}")
        cmd.create('hydrophilic_surf', f"{obj_name} and resn {'+'.join(hydrophilic)}")
        cmd.show('surface', 'hydrophobic_surf')
        cmd.show('surface', 'hydrophilic_surf')
        cmd.color('hydrophobic_color', 'hydrophobic_surf')
        cmd.color('hydrophilic_color', 'hydrophilic_surf')
        cmd.set('transparency', 0.5, 'hydrophobic_surf')
        cmd.set('transparency', 0.5, 'hydrophilic_surf')
    elif mode == 'highlight':
        cmd.show('sticks', f"{obj_name} and resn {'+'.join(hydrophobic)} and sidechain")
        cmd.show('sticks', f"{obj_name} and resn {'+'.join(hydrophilic)} and sidechain")
        cmd.color('hydrophobic_color', f"{obj_name} and resn {'+'.join(hydrophobic)} and sidechain")
        cmd.color('hydrophilic_color', f"{obj_name} and resn {'+'.join(hydrophilic)} and sidechain")

def find_glycosylation_sites(obj_name):
    """
    Find N-glycosylation sites (N-X-S/T pattern) and actual glycosylation
    
    Args:
        obj_name: PyMOL object name
    
    Returns:
        List of (chain, resi, type) tuples where type is 'potential' or 'actual'
    """
    glyco_sites = []
    
    # Get all residues
    stored.residues = []
    cmd.iterate(f"{obj_name} and name CA", "stored.residues.append((chain, resi, resn))")
    
    # Group by chain
    chains = {}
    for chain, resi, resn in stored.residues:
        if chain not in chains:
            chains[chain] = []
        chains[chain].append((int(resi), resn))
    
    # Find N-X-S/T pattern
    for chain, residues in chains.items():
        residues.sort()
        for i in range(len(residues) - 2):
            resi1, resn1 = residues[i]
            resi2, resn2 = residues[i + 1]
            resi3, resn3 = residues[i + 2]
            
            if resn1 == 'ASN' and resn2 != 'PRO' and resn3 in ['SER', 'THR']:
                if resi2 == resi1 + 1 and resi3 == resi2 + 1:
                    glyco_sites.append((chain, str(resi1), 'potential'))
    
    # Check for actual glycosylation (presence of sugar residues)
    sugar_residues = ['NAG', 'MAN', 'FUC', 'GAL', 'BMA', 'NDG']
    stored.sugar_nearby = []
    for chain, resi, gtype in glyco_sites:
        nearby = cmd.select('temp_sugar', f"resn {'+'.join(sugar_residues)} within 5 of ({obj_name} and chain {chain} and resi {resi})")
        if nearby > 0:
            glyco_sites[glyco_sites.index((chain, resi, gtype))] = (chain, resi, 'actual')
        cmd.delete('temp_sugar')
    
    return glyco_sites

def mark_glycosylation_sites(obj_name, sites, mode='color'):
    """
    Mark glycosylation sites
    
    Args:
        obj_name: PyMOL object name
        sites: List of (chain, resi, type) tuples
        mode: 'color', 'surface', or 'highlight'
    """
    if not sites:
        return
    
    # Define custom colors
    cmd.set_color('glyco_potential_color', [1.0, 0.4, 0.8])  # bright pink
    cmd.set_color('glyco_actual_color', [0.8, 0.0, 0.5])     # deep pink
    
    for chain, resi, gtype in sites:
        sel = f"{obj_name} and chain {chain} and resi {resi}"
        color = 'glyco_actual_color' if gtype == 'actual' else 'glyco_potential_color'
        
        if mode == 'color':
            cmd.color(color, sel)
        elif mode == 'surface':
            cmd.create(f'glyco_{chain}_{resi}', sel)
            cmd.show('surface', f'glyco_{chain}_{resi}')
            cmd.color(color, f'glyco_{chain}_{resi}')
            cmd.set('transparency', 0.5, f'glyco_{chain}_{resi}')
        elif mode == 'highlight':
            cmd.show('sticks', sel)
            cmd.color(color, sel)
            cmd.set('stick_radius', 0.3, sel)

def visualize_molecule(pdb_file, output_dir=None, show_sidechain=False, quality='low', view_style='cartoon', 
                       analysis_mode='cysteine'):
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
    
    # Color chains based on analysis mode
    if analysis_mode == 'cysteine':
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
    
    elif analysis_mode == 'property':
        # Use light gray as base color for property analysis
        cmd.color("gray70", "molecule")
    
    # Apply view style
    cmd.show(view_style, "molecule")
    
    # Analysis mode: cysteine or property
    free_cys = []
    glyco_sites = []
    
    if analysis_mode == 'cysteine':
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
    
    elif analysis_mode == 'property':
        # Mark hydrophobic/hydrophilic residues
        mark_hydrophobic_hydrophilic("molecule", 'color')
        
        # Find and mark glycosylation sites
        glyco_sites = find_glycosylation_sites("molecule")
        # Use color mode for surface view, highlight mode for others
        glyco_mode = 'color' if view_style == 'surface' else 'highlight'
        mark_glycosylation_sites("molecule", glyco_sites, glyco_mode)
    
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
    
    # Generate close-ups based on analysis mode
    if analysis_mode == 'cysteine' and free_cys:
        # Free cysteine close-ups
        for idx, (chain, resi, min_dist) in enumerate(free_cys):
            cmd.orient(f"molecule and chain {chain} and resi {resi}")
            cmd.zoom(f"molecule and chain {chain} and resi {resi}", 15)
            cmd.label(f"molecule and chain {chain} and resi {resi} and name CA", f"'Chain {chain} CYS {resi}'")
            cmd.set("label_size", 20)
            cmd.set("label_color", "red")
            cmd.png(f"{output_dir}/{base_name}_free_cys_{idx+1}.png", width=1920, height=1080, dpi=300, ray=1)
            cmd.hide("labels")
    
    elif analysis_mode == 'property' and glyco_sites:
        # Glycosylation site close-ups
        for idx, (chain, resi, gtype) in enumerate(glyco_sites):
            cmd.orient(f"molecule and chain {chain} and resi {resi}")
            cmd.zoom(f"molecule and chain {chain} and resi {resi}", 15)
            label_text = f"'Chain {chain} ASN {resi} ({gtype})'"
            cmd.label(f"molecule and chain {chain} and resi {resi} and name CA", label_text)
            cmd.set("label_size", 20)
            cmd.set("label_color", "glyco_actual_color" if gtype == 'actual' else "glyco_potential_color")
            cmd.png(f"{output_dir}/{base_name}_glyco_site_{idx+1}.png", width=1920, height=1080, dpi=300, ray=1)
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
    
    # Zoom to target based on analysis mode
    zoom_target = None
    if analysis_mode == 'cysteine' and free_cys:
        chain, resi, min_dist = free_cys[0]
        zoom_target = (chain, resi)
    elif analysis_mode == 'property' and glyco_sites:
        chain, resi, gtype = glyco_sites[0]
        zoom_target = (chain, resi)
    
    if zoom_target:
        chain, resi = zoom_target
        
        # Save the current view after rotation
        view_start = cmd.get_view()
        
        # Get the target view (zoomed to target residue)
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
    
    if analysis_mode == 'cysteine':
        print(f"Free cysteines found: {len(free_cys)}")
        if free_cys:
            for chain, resi, min_dist in free_cys:
                dist_str = f"{min_dist:.2f} Å" if min_dist else "N/A"
                print(f"  - Chain {chain}, Residue {resi} (nearest SG: {dist_str})")
    
    elif analysis_mode == 'property':
        print(f"Glycosylation sites found: {len(glyco_sites)}")
        if glyco_sites:
            for chain, resi, gtype in glyco_sites:
                print(f"  - Chain {chain}, Asn {resi} ({gtype})")
    
    print(f"\nTo view interactively: pymol {session_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Visualize protein structures with PyMOL with two analysis modes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Description:
  This script loads a PDB file and visualizes the protein structure using PyMOL.
  It supports two analysis modes:
  
  1. Cysteine Analysis (default): Identifies free cysteines not involved in disulfide bonds
  2. Property Analysis: Marks hydrophobic/hydrophilic residues and N-glycosylation sites

Analysis Modes:
  cysteine  - Free cysteine detection for ADC development
              • Marks all cysteine residues (YELLOW sticks)
              • Highlights free cysteines (RED sticks, thicker)
              • Generates close-up views of each free cysteine
              • GIF animation zooms to first free cysteine
              • Chains colored with Mirabo palette (multi-color)
              • Reports S-S distances for each free cysteine
  
  property  - Molecular property analysis for surface characterization
              • Marks hydrophobic residues (tan/ochre color)
              • Marks hydrophilic residues (steel blue color)
              • Highlights N-glycosylation sites (pink, display adapts to view mode)
                - Surface view: colored surface (shows surface accessibility)
                - Other views: sticks (shows precise residue position)
              • Generates close-up views of each glycosylation site
              • GIF animation zooms to first glycosylation site
              • Chains colored in light gray (unified base color)
              • Detects N-X-S/T pattern and actual glycosylation

View Styles:
  cartoon  - Secondary structure representation showing alpha helices and beta sheets
             Best for: Understanding overall protein fold and domain architecture
  
  surface  - Molecular surface representation showing the solvent-accessible surface
             Best for: Analyzing binding sites, protein-protein interfaces, cavities
             Recommended with: --analysis-mode property (shows hydrophobic distribution)
  
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

Color Schemes:
  Cysteine Analysis Mode:
    • Chains: Mirabo-inspired palette (deep blue, lime green, orange, purple, cyan, gold, magenta, light blue)
    • All cysteines: YELLOW sticks
    • Free cysteines: RED sticks (thicker, potential ADC conjugation sites)
  
  Property Analysis Mode:
    • Chains: Light gray (unified base color to reduce visual clutter)
    • Hydrophobic residues (ALA,VAL,LEU,ILE,PHE,TRP,MET,PRO): Tan/ochre
    • Hydrophilic residues (SER,THR,ASN,GLN,TYR,CYS,LYS,ARG,HIS,ASP,GLU): Steel blue
    • Glycosylation sites (display adapts to view mode):
      - Surface view: Colored surface (potential: bright pink, actual: deep pink)
      - Other views: Sticks (potential: bright pink, actual: deep pink)

Output Files:
  <name>_front.png       - Front view of the molecule
  <name>_top.png         - Top view (90° rotation on X-axis)
  <name>_side.png        - Side view (90° rotation on Y-axis)
  <name>_overview.png    - Zoomed out overview
  
  Cysteine mode:
    <name>_free_cys_N.png  - Close-up of each free cysteine (if found)
    <name>_rotation.gif    - 360° rotation + zoom to first free cysteine
  
  Property mode:
    <name>_glyco_site_N.png - Close-up of each glycosylation site (if found)
    <name>_rotation.gif     - 360° rotation + zoom to first glycosylation site
  
  <name>_session.pse     - PyMOL session file for interactive exploration

Examples:
  # Basic usage - cysteine analysis (default)
  %(prog)s protein.pdb
  
  # Property analysis - hydrophobic/hydrophilic + glycosylation
  %(prog)s protein.pdb --analysis-mode property
  
  # Surface view for property analysis (recommended)
  %(prog)s protein.pdb --analysis-mode property --view surface
  
  # High quality output
  %(prog)s protein.pdb output/ --analysis-mode property --view surface --quality high
  
  # Cysteine analysis with surface view
  %(prog)s protein.pdb output/ --view surface --quality high
  
  # Show side chains for detailed analysis
  %(prog)s protein.pdb --show-sidechain --view cartoon

Tips:
  - Use 'cysteine' mode for ADC development and disulfide bond analysis
  - Use 'property' mode for glycosylation engineering and surface characterization
  - Combine 'property' mode with 'surface' view for best hydrophobic visualization
  - Use 'cartoon' for presentations and overall structure analysis
  - Open the .pse session file in PyMOL for interactive manipulation
  - For IgG antibodies, property mode will detect Asn297 glycosylation site in Fc region
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
    
    parser.add_argument('--analysis-mode', 
                        choices=['cysteine', 'property'],
                        default='cysteine',
                        metavar='MODE',
                        help='Analysis mode: cysteine (default, mark free cysteines) or property (mark hydrophobic/hydrophilic and glycosylation sites)')
    
    args = parser.parse_args()
    
    pymol.finish_launching(['pymol', '-c'])
    visualize_molecule(args.pdb_file, args.output_dir, args.show_sidechain, args.quality, args.view, args.analysis_mode)
