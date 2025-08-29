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

    with open(pdb_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                is_in_atom_section = True
                chain_id = line[21]
                if chain_id not in atom_data_by_chain:
                    atom_data_by_chain[chain_id] = []
                atom_data_by_chain[chain_id].append(line)
            elif is_in_atom_section:
                # Once we've seen atoms, everything else is a footer until END
                if not line.strip().startswith('END'):
                     footer_lines.append(line)
            else:
                header_lines.append(line)

    # --- Pass 2: Determine the processing order and create a renumbering plan ---
    processing_order = []
    processed_chains = set()
    if chain_groups:
        for group in chain_groups:
            for chain_id in group:
                if chain_id in atom_data_by_chain:
                    processing_order.append(chain_id)
                    processed_chains.add(chain_id)
    # Add any chains from the PDB that weren't specified in the info file
    for chain_id in atom_data_by_chain:
        if chain_id not in processed_chains:
            processing_order.append(chain_id)

    # --- Pass 3: Renumber and write the new file in the correct order ---
    chain_info = OrderedDict()
    residue_counter = 0
    
    with open(output_pdb_file, 'w') as f_out:
        # Write the original file's header
        f_out.writelines(header_lines)

        # Write the ATOM/HETATM data in the specified order
        for chain_id in processing_order:
            if chain_id in atom_data_by_chain:
                last_original_res_num = None
                
                for line in atom_data_by_chain[chain_id]:
                    # Only renumber ATOM records, leave HETATM as is
                    if line.startswith('ATOM'):
                        original_res_num = line[22:26].strip()
                        
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
                
                # Add a TER card to signify the end of the chain's ATOM records
                f_out.write("TER\n")
        
        # Write the original file's footer (CONECT, MASTER, etc.)
        f_out.writelines(footer_lines)
        f_out.write("END\n")

    return chain_info

def parse_info_file(info_file):
    """Parses the information file to get the chain groupings."""
    chain_groups = []
    with open(info_file, 'r') as f:
        for line in f:
            if 'protein' in line:
                parts = line.split(':')
                if len(parts) > 1:
                    chains = parts[1].strip().split(',')
                    chain_groups.append([chain.strip() for chain in chains])
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
    else:
        print("No 'group_information.txt' found. Renumbering chains sequentially.")

    temp_output_pdb_path = pdb_file_path + ".tmp"

    chain_info = renumber_pdb_structure_preserved(pdb_file_path, temp_output_pdb_path, chain_groups)

    os.replace(temp_output_pdb_path, pdb_file_path)

    info_dir = "residual_number_information"
    os.makedirs(info_dir, exist_ok=True)
    
    summary_file_path = os.path.join(info_dir, f"{folder_name}_{base_name}.txt")

    with open(summary_file_path, 'w') as f_summary:
        f_summary.write(f"Summary for {pdb_file_name}:\n")
        for chain_id, info in chain_info.items():
            f_summary.write(f"Chain {chain_id}: {info['start']}-{info['end']}\n")

    print(f"Successfully renumbered and replaced: {pdb_file_path}")
    print(f"Chain summary file created: {summary_file_path}")

if __name__ == "__main__":
    main()
