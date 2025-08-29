import os
import sys
from collections import OrderedDict

def renumber_pdb_structure_preserved(pdb_file, output_pdb_file, chain_groups=None):
    """
    Renumbers residues in a PDB file continuously while perfectly preserving
    the original file structure, including HETATM and other records.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_pdb_file (str): Path to the temporary output PDB file.
        chain_groups (list of list of str, optional): Defines the order for
            continuous numbering. Defaults to None.

    Returns:
        OrderedDict: A dictionary with the start and end residue numbers for
                     each chain after renumbering.
    """
    # --- Pass 1: Read and categorize every line of the PDB file ---
    header_lines = []
    footer_lines = []
    # Store ATOM and HETATM lines together to maintain their relative order
    atom_data_by_chain = OrderedDict()
    is_in_atom_section = False

    try:
        with open(pdb_file, 'r') as f_in:
            line_num = 0
            for line in f_in:
                line_num += 1
                try:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Check if line is long enough to contain chain ID
                        if len(line) < 22:
                            print(f"Warning: Line {line_num} is too short to contain chain ID: {line.strip()}")
                            continue
                            
                        is_in_atom_section = True
                        chain_id = line[21]
                        
                        # Handle empty or space chain IDs
                        if chain_id == ' ':
                            chain_id = 'A'  # Default to chain A for unnamed chains
                            print(f"Warning: Empty chain ID on line {line_num}, defaulting to chain 'A'")
                            
                        if chain_id not in atom_data_by_chain:
                            atom_data_by_chain[chain_id] = []
                        atom_data_by_chain[chain_id].append(line)
                    elif is_in_atom_section:
                        # Once we've seen atoms, everything else is a footer until END
                        if not line.strip().startswith('END'):
                             footer_lines.append(line)
                    else:
                        header_lines.append(line)
                except (IndexError, UnicodeDecodeError) as e:
                    print(f"Warning: Error processing line {line_num} in {pdb_file}: {e}")
                    print(f"  Problematic line: {line.strip()}")
                    continue
                    
    except FileNotFoundError:
        print(f"Error: PDB file not found: {pdb_file}")
        return OrderedDict()
    except PermissionError:
        print(f"Error: Permission denied reading file: {pdb_file}")
        return OrderedDict()
    except UnicodeDecodeError:
        print(f"Error: Unable to decode file (encoding issue): {pdb_file}")
        return OrderedDict()
    except Exception as e:
        print(f"Error: Unexpected error reading {pdb_file}: {e}")
        return OrderedDict()

    # --- Pass 2: Determine the processing order and create a renumbering plan ---
    processing_order = []
    processed_chains = set()
    if chain_groups:
        for group in chain_groups:
            for chain_id in group:
                if chain_id in atom_data_by_chain:
                    processing_order.append(chain_id)
                    processed_chains.add(chain_id)
                else:
                    print(f"Warning: Chain '{chain_id}' specified in group info but not found in PDB file")
    # Add any chains from the PDB that weren't specified in the info file
    for chain_id in atom_data_by_chain:
        if chain_id not in processed_chains:
            processing_order.append(chain_id)

    # Check if we have any chains to process
    if not processing_order:
        print(f"Error: No valid chains found in PDB file {pdb_file}")
        return OrderedDict()

    # --- Pass 3: Renumber and write the new file in the correct order ---
    chain_info = OrderedDict()
    residue_counter = 0
    
    try:
        with open(output_pdb_file, 'w') as f_out:
            # Write the original file's header
            try:
                f_out.writelines(header_lines)
            except Exception as e:
                print(f"Warning: Error writing header lines: {e}")

            # Write the ATOM/HETATM data in the specified order
            for chain_id in processing_order:
                if chain_id in atom_data_by_chain:
                    last_original_res_num = None
                    
                    for line in atom_data_by_chain[chain_id]:
                        try:
                            # Only renumber ATOM records, leave HETATM as is
                            if line.startswith('ATOM'):
                                # Check if line is long enough to contain residue number
                                if len(line) < 26:
                                    print(f"Warning: ATOM line too short to contain residue number: {line.strip()}")
                                    f_out.write(line)  # Write as-is
                                    continue
                                    
                                original_res_num = line[22:26].strip()
                                
                                if not original_res_num:
                                    print(f"Warning: Empty residue number in line: {line.strip()}")
                                    f_out.write(line)  # Write as-is
                                    continue
                                
                                if original_res_num != last_original_res_num:
                                    residue_counter += 1
                                    last_original_res_num = original_res_num
                                
                                if chain_id not in chain_info:
                                    chain_info[chain_id] = {'start': residue_counter}
                                
                                chain_info[chain_id]['end'] = residue_counter
                                
                                new_line = line[:22] + str(residue_counter).rjust(4) + line[26:]
                                f_out.write(new_line)
                            else:
                                # Write HETATM lines without modification
                                f_out.write(line)
                        except Exception as e:
                            print(f"Warning: Error processing line for chain {chain_id}: {e}")
                            print(f"  Problematic line: {line.strip()}")
                            f_out.write(line)  # Write original line as fallback
                            continue
                    
                    # Add a TER card to signify the end of the chain's ATOM records
                    try:
                        f_out.write("TER\n")
                    except Exception as e:
                        print(f"Warning: Error writing TER record for chain {chain_id}: {e}")
            
            # Write the original file's footer (CONECT, MASTER, etc.)
            try:
                f_out.writelines(footer_lines)
                f_out.write("END\n")
            except Exception as e:
                print(f"Warning: Error writing footer lines: {e}")
                
    except PermissionError:
        print(f"Error: Permission denied writing to file: {output_pdb_file}")
        return OrderedDict()
    except IOError as e:
        print(f"Error: I/O error writing to file {output_pdb_file}: {e}")
        return OrderedDict()
    except Exception as e:
        print(f"Error: Unexpected error writing to file {output_pdb_file}: {e}")
        return OrderedDict()

    return chain_info

def parse_info_file(info_file):
    """Parses the information file to get the chain groupings."""
    chain_groups = []
    try:
        with open(info_file, 'r') as f:
            line_num = 0
            for line in f:
                line_num += 1
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                    
                if ':' in line:
                    try:
                        parts = line.split(':')
                        if len(parts) > 1:
                            group_name = parts[0].strip()
                            chain_list = parts[1].strip()
                            
                            if not group_name:
                                print(f"Warning: Empty group name on line {line_num} in {info_file}")
                                continue
                                
                            if not chain_list:
                                print(f"Warning: Empty chain list for group '{group_name}' on line {line_num}")
                                continue
                                
                            chains = [chain.strip() for chain in chain_list.split(',') if chain.strip()]
                            if chains:
                                chain_groups.append(chains)
                            else:
                                print(f"Warning: No valid chains found for group '{group_name}' on line {line_num}")
                        else:
                            print(f"Warning: Malformed line (no content after ':') on line {line_num} in {info_file}: {line}")
                    except Exception as e:
                        print(f"Warning: Error parsing line {line_num} in {info_file}: {line}")
                        print(f"  Error details: {e}")
                        continue
                else:
                    # Line doesn't contain ':', might be malformed
                    if line:  # Only warn for non-empty lines
                        print(f"Warning: Skipping malformed line {line_num} in {info_file}: {line}")
                        
    except FileNotFoundError:
        print(f"Error: Group information file not found: {info_file}")
        return []
    except PermissionError:
        print(f"Error: Permission denied reading file: {info_file}")
        return []
    except UnicodeDecodeError:
        print(f"Error: Unable to decode file (encoding issue): {info_file}")
        return []
    except Exception as e:
        print(f"Error: Unexpected error reading {info_file}: {e}")
        return []
        
    return chain_groups

def main():
    """Main function to run the script from the command line."""
    if len(sys.argv) != 2:
        print("Usage: python renumber_residual.py <folder_name>")
        sys.exit(1)

    folder_name = sys.argv[1]

    if not os.path.isdir(folder_name):
        print(f"Error: Folder '{folder_name}' not found in the current directory.")
        sys.exit(1)

    # --- Start of Modified Section ---

    # Construct the specific filename based on the folder name
    # os.path.basename handles cases like './folder' or 'folder/'
    #base_folder_name = os.path.basename(os.path.normpath(folder_name))
    pdb_file_name = f"solvated.pdb"
    pdb_file_path = os.path.join(folder_name, pdb_file_name)

    # Check if the specific file exists. If not, exit.
    if not os.path.isfile(pdb_file_path):
        print(f"Error: Required file not found at '{pdb_file_path}'")
        sys.exit(1)

    # --- End of Modified Section ---

    print(f"--- Processing: {pdb_file_name} ---")
    base_name = os.path.splitext(pdb_file_name)[0]

    info_file_path = os.path.join(folder_name, 'group_information.txt')
    chain_groups = None
    if os.path.isfile(info_file_path):
        print(f"Found '{info_file_path}', using it for chain grouping.")
        chain_groups = parse_info_file(info_file_path)
        
        # Check if parsing was successful
        if not chain_groups:
            print("Warning: No valid groups found in group information file. Proceeding with sequential chain ordering.")
            chain_groups = None
    else:
        print("No 'group_information.txt' found. Renumbering chains sequentially.")

    temp_output_pdb_path = pdb_file_path + ".tmp"

    # Attempt to renumber the PDB file
    try:
        chain_info = renumber_pdb_structure_preserved(pdb_file_path, temp_output_pdb_path, chain_groups)
        
        # Check if renumbering was successful
        if not chain_info:
            print(f"Error: Failed to renumber PDB file '{pdb_file_path}'. No output generated.")
            # Clean up temporary file if it exists
            if os.path.isfile(temp_output_pdb_path):
                try:
                    os.remove(temp_output_pdb_path)
                except Exception as e:
                    print(f"Warning: Could not remove temporary file {temp_output_pdb_path}: {e}")
            sys.exit(1)
            
        # Replace original file with renumbered version
        try:
            os.replace(temp_output_pdb_path, pdb_file_path)
        except Exception as e:
            print(f"Error: Failed to replace original file with renumbered version: {e}")
            # Clean up temporary file
            if os.path.isfile(temp_output_pdb_path):
                try:
                    os.remove(temp_output_pdb_path)
                except Exception as cleanup_e:
                    print(f"Warning: Could not remove temporary file {temp_output_pdb_path}: {cleanup_e}")
            sys.exit(1)

    except Exception as e:
        print(f"Error: Unexpected error during renumbering: {e}")
        # Clean up temporary file if it exists
        if os.path.isfile(temp_output_pdb_path):
            try:
                os.remove(temp_output_pdb_path)
            except Exception as cleanup_e:
                print(f"Warning: Could not remove temporary file {temp_output_pdb_path}: {cleanup_e}")
        sys.exit(1)

    # Create summary file
    info_dir = "residual_number_information"
    try:
        os.makedirs(info_dir, exist_ok=True)
    except Exception as e:
        print(f"Error: Could not create directory '{info_dir}': {e}")
        sys.exit(1)
    
    summary_file_path = os.path.join(info_dir, f"{folder_name}_{base_name}.txt")

    try:
        with open(summary_file_path, 'w') as f_summary:
            f_summary.write(f"Summary for {pdb_file_name}:\n")
            for chain_id, info in chain_info.items():
                f_summary.write(f"Chain {chain_id}: {info['start']}-{info['end']}\n")
    except Exception as e:
        print(f"Error: Could not write summary file '{summary_file_path}': {e}")
        sys.exit(1)

    print(f"Successfully renumbered and replaced: {pdb_file_path}")
    print(f"Chain summary file created: {summary_file_path}")
    print(f"Total chains processed: {len(chain_info)}")
    print("Chain information:")
    for chain_id, info in chain_info.items():
        residue_count = info['end'] - info['start'] + 1
        print(f"  Chain {chain_id}: residues {info['start']}-{info['end']} ({residue_count} residues)")

if __name__ == "__main__":
    main()
