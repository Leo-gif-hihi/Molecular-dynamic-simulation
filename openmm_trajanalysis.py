#!/usr/bin/env python3
"""
openmm_trajanalysis.py

Compute and plot:
  0) Detected total simulation length
  1) Backbone RMSD vs. time
  2) Radius of gyration vs. time
  3) Cα RMSF per residue

This script will automatically perform a grouped RMSF analysis if it finds both:
  - <pdbid>/group_information.txt
  - residual_number_information/<pdbid>_solvated.txt
If these files are found, it will generate a single plot with each group colored
differently. By default, it saves all grouped data to a single text file.
Use the --split-rmsf flag to save a separate data file for each group instead.

Assumes you have:
  <pdbid>/solvated.pdb
  <pdbid>/prod_full.dcd

Usage:
    python3 openmm_trajanalysis.py <pdbid> [options]

Positional arguments:
  pdbid                    4-letter PDB ID directory (e.g. 3e7y)

Optional arguments:
  -t, --topology PATH      Topology PDB (default: <pdbid>/solvated.pdb)
  -x, --trajectory PATH    Trajectory DCD (default: <pdbid>/prod_full.dcd)
  -i, --interval FLOAT     Frame interval in ps (default: detect from DCD header)
  -o, --outdir DIR         Output directory for plots
  --split-rmsf             If group info files are found, split RMSF data into
                           separate files for each protein group.
"""

import os
os.environ['MPLBACKEND'] = 'Agg'  # <-- Set the backend BEFORE any other imports

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis import rms
import mdtraj as md

def detect_interval_ps(u):
    """Try to read u.trajectory.ts.dt (ps); fallback to None."""
    try:
        return float(u.trajectory.ts.dt)
    except Exception:
        return None

def parse_group_info(filepath):
    """Parses group_information.txt into a dictionary like {'protein 1': ['A', 'C']}. """
    groups = {}
    try:
        with open(filepath, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                    
                if ':' in line:
                    try:
                        key, val = line.split(':', 1)
                        key = key.strip()
                        
                        if not key:
                            print(f"Warning: Empty group name on line {line_num} in {filepath}")
                            continue
                            
                        # Handle empty or whitespace-only values
                        val = val.strip()
                        if not val:
                            print(f"Warning: Empty chain list for group '{key}' on line {line_num}")
                            groups[key] = []
                            continue
                            
                        # Parse chains and filter out empty ones
                        chains = [c.strip() for c in val.split(',') if c.strip()]
                        if not chains:
                            print(f"Warning: No valid chains found for group '{key}' on line {line_num}")
                            groups[key] = []
                        else:
                            groups[key] = chains
                            
                    except ValueError as e:
                        print(f"Warning: Error parsing line {line_num} in {filepath}: {line}")
                        print(f"  Error details: {e}")
                        continue
                else:
                    # Line doesn't contain ':', might be malformed
                    if line:  # Only warn for non-empty lines
                        print(f"Warning: Skipping malformed line {line_num} in {filepath}: {line}")
                        
    except FileNotFoundError:
        print(f"Error: Group information file not found: {filepath}")
        return {}
    except PermissionError:
        print(f"Error: Permission denied reading file: {filepath}")
        return {}
    except UnicodeDecodeError:
        print(f"Error: Unable to decode file (encoding issue): {filepath}")
        return {}
    except Exception as e:
        print(f"Error: Unexpected error reading {filepath}: {e}")
        return {}
        
    return groups

def parse_residue_info(filepath):
    """Parses <pdbid>_solvated.txt for chain-to-residue ranges into a dict like {'A': (1, 21)}. """
    chain_ranges = {}
    try:
        with open(filepath, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                    
                if line.startswith('Chain'):
                    try:
                        parts = line.split()
                        if len(parts) < 3:
                            print(f"Warning: Incomplete chain information on line {line_num} in {filepath}: {line}")
                            continue
                            
                        chain_id = parts[1].replace(':', '')
                        if not chain_id:
                            print(f"Warning: Empty chain ID on line {line_num} in {filepath}")
                            continue
                            
                        range_str = parts[2]
                        if '-' not in range_str:
                            print(f"Warning: Invalid range format on line {line_num} in {filepath}: {range_str}")
                            continue
                            
                        try:
                            start, end = map(int, range_str.split('-'))
                            if start > end:
                                print(f"Warning: Invalid range (start > end) on line {line_num} in {filepath}: {start}-{end}")
                                continue
                            chain_ranges[chain_id] = (start, end)
                        except ValueError as e:
                            print(f"Warning: Could not parse residue range on line {line_num} in {filepath}: {range_str}")
                            print(f"  Error details: {e}")
                            continue
                            
                    except (IndexError, ValueError) as e:
                        print(f"Warning: Error parsing chain line {line_num} in {filepath}: {line}")
                        print(f"  Error details: {e}")
                        continue
                        
    except FileNotFoundError:
        print(f"Error: Residue information file not found: {filepath}")
        return {}
    except PermissionError:
        print(f"Error: Permission denied reading file: {filepath}")
        return {}
    except UnicodeDecodeError:
        print(f"Error: Unable to decode file (encoding issue): {filepath}")
        return {}
    except Exception as e:
        print(f"Error: Unexpected error reading {filepath}: {e}")
        return {}
        
    return chain_ranges

def parse_args():
    p = argparse.ArgumentParser(description="Analyze merged trajectory for a pdbid")
    p.add_argument("pdbid", help="4-letter PDB ID directory")
    p.add_argument("-t", "--topology",
                   help="Topology PDB (default: <pdbid>/solvated.pdb)")
    p.add_argument("-x", "--trajectory",
                   help="Trajectory DCD (default: <pdbid>/prod_full.dcd)")
    p.add_argument("-i", "--interval", type=float,
                   help="Frame interval in ps (default: detect from header)")
    p.add_argument("-o", "--outdir",
                   help="Directory to save plots (default: analysis_<pdbid>_<timestamp>)")
    p.add_argument("--split-rmsf", action="store_true",
                   help="If group info files are found, split RMSF data into "
                        "separate files for each protein group.")
    return p.parse_args()

def main():
    args = parse_args()
    pdbid = args.pdbid

    # Resolve paths
    top_def   = os.path.join(pdbid, "solvated.pdb")
    traj_def  = os.path.join(pdbid, "prod_full.dcd")
    topo_path = os.path.abspath(args.topology or top_def)
    traj_path = os.path.abspath(args.trajectory or traj_def)

    if not os.path.isfile(topo_path):
        sys.exit(f"Topology file not found: {topo_path}")
    if not os.path.isfile(traj_path):
        sys.exit(f"Trajectory file not found: {traj_path}")

    # Create output directory
    outdir = args.outdir or pdbid
    os.makedirs(outdir, exist_ok=True)

    # Load MDAnalysis Universe
    print("Loading trajectory with MDAnalysis...")
    u = mda.Universe(topo_path, traj_path)
    n_frames = len(u.trajectory)
    # Determine frame spacing
    interval = args.interval or detect_interval_ps(u)
    if interval is None:
        sys.exit("Failed to detect frame interval; please specify --interval")
    total_ns = (n_frames - 1) * interval / 1000.0
    print(f"Detected {n_frames} frames, interval = {interval:.2f} ps → total ≈ {total_ns:.3f} ns\n")

    # Define time axis in nanoseconds
    times_ns = (np.arange(n_frames) * interval) / 1000.0

    # 1) RMSD (backbone)
    print("Computing RMSD...")
    rmsd_calc = rms.RMSD(u, u, select="backbone", ref_frame=0)
    rmsd_calc.run()
    rmsd_data = rmsd_calc.results.rmsd
    rmsd_nm   = rmsd_data[:,2] / 10.0

    plt.figure()
    plt.plot(times_ns, rmsd_nm, "-o", markersize=3)
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (nm)")
    plt.title("Backbone RMSD vs. Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rmsd_vs_time.png"), dpi=300)
    plt.close()
    print(f"  → Saved {outdir}/rmsd_vs_time.png")
    rmsd_output_path = os.path.join(outdir, "rmsd_data.txt")
    np.savetxt(rmsd_output_path,
               np.column_stack((times_ns, rmsd_nm)),
               header="Time(ns) RMSD(nm)",
               fmt="%.4f")
    print(f"  → Saved {rmsd_output_path}")
    
    # 2) Radius of gyration
    print("\nComputing Radius of Gyration...")
    heavy   = u.select_atoms("not name H*")
    rg_vals = []
    for ts in u.trajectory:
        rg_vals.append(heavy.radius_of_gyration() / 10.0)
    rg_vals = np.array(rg_vals)

    plt.figure()
    plt.plot(times_ns, rg_vals, "-o", markersize=3)
    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (nm)")
    plt.title("Rg vs. Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "rg_vs_time.png"), dpi=300)
    plt.close()
    print(f"  → Saved {outdir}/rg_vs_time.png")
    rg_output_path = os.path.join(outdir, "rg_data.txt")
    np.savetxt(rg_output_path,
               np.column_stack((times_ns, rg_vals)),
               header="Time(ns) Rg(nm)",
               fmt="%.4f")
    print(f"  → Saved {rg_output_path}")
    
    # 3) RMSF (MDTraj)
    print("\nComputing RMSF...")
    traj = md.load(traj_path, top=topo_path)
    traj.superpose(traj, 0)
    ca_idx  = traj.topology.select("name CA")
    rmsf_A  = md.rmsf(traj, traj, atom_indices=ca_idx)
    rmsf_nm = rmsf_A / 10.0
    resids  = np.array([traj.topology.atom(i).residue.resSeq for i in ca_idx])

    # --- UPDATED RMSF ANALYSIS LOGIC ---
    group_info_path = os.path.join(pdbid, "group_information.txt")
    residue_info_path = os.path.join("residual_number_information", f"{pdbid}_solvated.txt")
    
    # If both info files are found, run a grouped analysis
    if os.path.isfile(group_info_path) and os.path.isfile(residue_info_path):
        print("  → Found info files. Performing grouped RMSF analysis.")
        chain_groups = parse_group_info(group_info_path)
        chain_ranges = parse_residue_info(residue_info_path)
        
        # Check if parsing was successful
        if not chain_groups:
            print("  → Warning: No valid groups found in group information file. Falling back to default analysis.")
            # Fall back to default non-grouped analysis
            rmsf_output_path = os.path.join(outdir, "rmsf_data.txt")
            np.savetxt(rmsf_output_path,
                       np.column_stack((resids, rmsf_nm)),
                       header="Residue_ID RMSF(nm)",
                       fmt=["%d", "%.4f"])
            print(f"  → Saved {rmsf_output_path}")
            
            plt.figure()
            plt.plot(resids, rmsf_nm, "-o", markersize=3)
            plt.xlabel("Residue")
            plt.ylabel("RMSF (nm)")
            plt.title("Backbone Cα RMSF per Residue")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"), dpi=300)
            plt.close()
            print(f"  → Saved {outdir}/rmsf_per_residue.png")
            print(f"\nAll plots saved in: {outdir}")
            return
            
        if not chain_ranges:
            print("  → Warning: No valid chain ranges found in residue information file. Falling back to default analysis.")
            # Fall back to default non-grouped analysis
            rmsf_output_path = os.path.join(outdir, "rmsf_data.txt")
            np.savetxt(rmsf_output_path,
                       np.column_stack((resids, rmsf_nm)),
                       header="Residue_ID RMSF(nm)",
                       fmt=["%d", "%.4f"])
            print(f"  → Saved {rmsf_output_path}")
            
            plt.figure()
            plt.plot(resids, rmsf_nm, "-o", markersize=3)
            plt.xlabel("Residue")
            plt.ylabel("RMSF (nm)")
            plt.title("Backbone Cα RMSF per Residue")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"), dpi=300)
            plt.close()
            print(f"  → Saved {outdir}/rmsf_per_residue.png")
            print(f"\nAll plots saved in: {outdir}")
            return

        plt.figure()
        
        # Check if the user wants to split the data files
        if args.split_rmsf:
            print("  → --split-rmsf flag detected. Splitting data output files by group.")
            for group_name, chains in chain_groups.items():
                group_residue_numbers = []
                for chain_id in chains:
                    if chain_id in chain_ranges:
                        start, end = chain_ranges[chain_id]
                        group_residue_numbers.extend(range(start, end + 1))
                    else:
                        print(f"Warning: Chain '{chain_id}' from group file not found in residue info file. Skipping.")

                mask = np.isin(resids, group_residue_numbers)
                current_resids = resids[mask]
                current_rmsf = rmsf_nm[mask]

                if len(current_resids) == 0:
                    print(f"Warning: No residues found for group '{group_name}'. Skipping.")
                    continue

                # Save a separate data file for this group
                sanitized_name = group_name.replace(" ", "_")
                output_path = os.path.join(outdir, f"rmsf_data_{sanitized_name}.txt")
                np.savetxt(output_path,
                           np.column_stack((current_resids, current_rmsf)),
                           header="Residue_ID RMSF(nm)",
                           fmt=["%d", "%.4f"])
                print(f"    → Saved {output_path}")

                plt.plot(current_resids, current_rmsf, "-o", markersize=3, label=group_name)
        
        else: # Default behavior: grouped plot, single combined data file
            print("  → No --split-rmsf flag. Generating single plot and single data file for all groups.")
            all_grouped_data = []
            for group_name, chains in chain_groups.items():
                group_residue_numbers = []
                for chain_id in chains:
                    if chain_id in chain_ranges:
                        start, end = chain_ranges[chain_id]
                        group_residue_numbers.extend(range(start, end + 1))
                    else:
                        print(f"Warning: Chain '{chain_id}' from group file not found in residue info file. Skipping.")

                mask = np.isin(resids, group_residue_numbers)
                current_resids = resids[mask]
                current_rmsf = rmsf_nm[mask]

                if len(current_resids) == 0:
                    print(f"Warning: No residues found for group '{group_name}'. Skipping.")
                    continue
                
                # Append data for combined text file
                group_names_col = np.full(current_resids.shape, group_name)
                all_grouped_data.append(np.column_stack((current_resids, current_rmsf, group_names_col)))

                plt.plot(current_resids, current_rmsf, "-o", markersize=3, label=group_name)
            
            # Save the combined data to a single file
            if all_grouped_data:
                combined_data = np.vstack(all_grouped_data)
                output_path = os.path.join(outdir, "rmsf_data_grouped.txt")
                np.savetxt(output_path,
                           combined_data,
                           header="Residue_ID RMSF(nm) Group",
                           fmt=['%s', '%.4f', '%s'])
                print(f"    → Saved combined data to {output_path}")

        # Finalize and save the plot for grouped analysis
        plt.xlabel("Residue")
        plt.ylabel("RMSF (nm)")
        plt.title("Backbone Cα RMSF per Residue (Grouped)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plot_path = os.path.join(outdir, "rmsf_per_residue_grouped.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"  → Saved grouped plot to {plot_path}")
    
    else: # Otherwise, run the default, non-grouped analysis
        print("  → Info file(s) for grouped analysis not found. Falling back to default.")
        if not os.path.isfile(group_info_path):
            print(f"    - Missing: {group_info_path}")
        if not os.path.isfile(residue_info_path):
            print(f"    - Missing: {residue_info_path}")
            
        rmsf_output_path = os.path.join(outdir, "rmsf_data.txt")
        np.savetxt(rmsf_output_path,
                   np.column_stack((resids, rmsf_nm)),
                   header="Residue_ID RMSF(nm)",
                   fmt=["%d", "%.4f"])
        print(f"  → Saved {rmsf_output_path}")
        
        plt.figure()
        plt.plot(resids, rmsf_nm, "-o", markersize=3)
        plt.xlabel("Residue")
        plt.ylabel("RMSF (nm)")
        plt.title("Backbone Cα RMSF per Residue")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "rmsf_per_residue.png"), dpi=300)
        plt.close()
        print(f"  → Saved {outdir}/rmsf_per_residue.png")

    print(f"\nAll plots saved in: {outdir}")

if __name__ == "__main__":
    main()