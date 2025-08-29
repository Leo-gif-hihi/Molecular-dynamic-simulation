-----

# Automated Molecular Dynamics Simulation Workflow for Protein Structures

This repository provides a streamlined workflow for running Molecular Dynamics (MD) simulations on protein structures using Google Colab. The process is automated through a series of Python scripts orchestrated by a Jupyter Notebook, leveraging powerful tools like OpenMM, PDBFixer, and MDAnalysis to prepare, simulate, and analyze protein PDB files. 🚀
## 📓 Colab Notebook 
You can view the **notebook** [here](https://colab.research.google.com/drive/16bE2P0ZposBOdU2zq9gWJJYCyOoFIEDn?usp=sharing).
## ✨ Features

  * **Automated Pipeline**: From PDB cleaning to final analysis, the entire workflow is handled by running a single notebook.
  * **Batch Processing**: Automatically processes multiple proteins, each located in its own subfolder.
  * **GPU Accelerated**: Utilizes Google Colab's free GPU resources for efficient OpenMM simulations.
  * **Flexible Analysis**: Generates key analysis plots for Root Mean Square Deviation (RMSD), Radius of Gyration (Rg), and Root Mean Square Fluctuation (RMSF).
  * **Group-Specific Analysis**: Supports optional, separate RMSF analysis for specific protein chains (e.g., for studying protein-protein interactions).
  * **Resumable Workflow**: Designed to be resumed from the last successfully processed protein if the Colab runtime is interrupted.

-----

## 📂 Project Structure

To use this workflow, your files must be organized in a specific structure within your Google Drive.

```
/your_main_folder/
│
├── openmm_proteinwater.py
├── openmm_trajanalysis.py
├── openmm_trajmerge.py
├── pdbfixer_cleaning.py
├── renumber_residual.py
│
├── /protein1_folder/
│   └── protein1.pdb
│
├── /protein2_folder/
│   └── protein2.pdb
│   └── group_information.txt   (Optional)
│
└── /protein3_folder/
    └── protein3.pdb

```

### **File Descriptions**

  * **Main Folder**: This is your main working directory (`workdir`). It must contain the five required Python scripts.
  * **Subfolders (`/protein1_folder/`, etc.)**:
      * Each subfolder represents a separate protein simulation.
      * It **must** contain exactly one PDB file (`.pdb`).
      * **Important**: The folder name must be all **lowercase**, with **no spaces** or **special characters**.
  * **`group_information.txt`** (Optional):
      * This file enables separate RMSF analysis for different groups of protein chains.
      * Place it inside any subfolder where you want this analysis performed.

-----

## 🚀 Setup and Usage

Follow these steps to set up and run the simulation workflow.

### 1\. Prepare Your Files in Google Drive

1.  **Create a main folder** in your Google Drive (e.g., `md_project`).
2.  **Upload the five Python scripts** (`pdbfixer_cleaning.py`, `openmm_proteinwater.py`, etc.) into this main folder. (have code perform this in the google colab notebook which ultilize the wget command)
3.  **Create subfolders** inside the main folder for each protein you want to simulate (e.g., `insulin`, `myoglobin`).
4.  **Place one PDB file** inside each of these subfolders.

### 2\. (Optional) Configure Group Analysis

If you want to perform separate RMSF analysis on specific chains (e.g., to analyze a receptor and a ligand separately), create a file named exactly `group_information.txt` inside the relevant protein subfolder and turn on --split-rmsf option (which have been done in the provided google colab notebook).

The format should be:

```
<name protein 1>: A,C
<name protein 2>: B,D
#Example:
insulin: A,C
alpha-glucosidase: B,D
```

This example will create two groups for RMSF analysis: one with chains A and C, and another with chains B and D. (view the example output in the output section)

### 3\. Configure and Run the Colab Notebook

1.  **Open the `md_simulation.ipynb` notebook** in Google Colab.
2.  **Mount your Google Drive** by running the first code cell.
3.  **Install Dependencies**: Run the second code cell to install Conda and all required packages (OpenMM, MDAnalysis, etc.). **The session will restart automatically after this step.**
4.  **Set the Working Directory**: In the third code cell, modify the `workdir` variable to point to the **main folder** you created in Google Drive.
    ```python
    # Define the path to your desired working folder in Google Drive
    workdir = '/content/drive/MyDrive/your_main_folder'
    ```
5.  **Run the Simulation**: Execute the final large code cell. The script will automatically find your protein subfolders and process them one by one. It can be safely re-run to resume processing if the runtime is disconnected.

-----

## 📊 Output

All main results are saved within a new folder named `analysis_{pdb_id}` created inside each protein's subfolder.

The main output includes:

  * **Analysis Data**: Text files (`.txt`) containing the raw data for RMSD, Rg, and RMSF.
  * **Plots**: PNG images visualizing the analysis results.
      * `rmsd_vs_time.png`
      * `rg_vs_time.png`
      * `rmsf_per_residue.png` (or `rmsf_per_residue_grouped.png` if using group analysis)
  * **Trajectory Files**: Merged trajectory (`prod_full.dcd`) and log (`prod_full.log`) files are stored in the protein's root subfolder.

You can view an example of a complete output folder [here](https://drive.google.com/drive/folders/1K_OSHJlJqwRNU5jvd-kDEdcdll42TbJS?usp=sharing).

-----

## 🔧 Troubleshooting

  * **Multiple PDB Files Error**:

      * **Error Message**: `Error: More than one PDB file found in 'folder_name'.`
      * **Solution**: Ensure only **one** `.pdb` file is present in each protein subfolder before starting the simulation. The script generates intermediate files like `_cleaned.pdb` and `solvated.pdb`, but the initial check is for a single input PDB.

  * **OpenMM Checkpoint Error**:

      * **Error Message**: `openmm.OpenMMException: loadCheckpoint: Checkpoint header was not correct`
      * **Solution**: This can happen if a simulation is interrupted and the checkpoint file (`prod.chk`) becomes corrupted. To fix this, delete the `prod.chk` file from the folder of the protein that failed and re-run the main script. The simulation for that protein will restart.

-----

## 📜 Citation

> Paul, S. K., Saddam, M., Tabassum, N., & Hasan, M. (2024). Molecular dynamics simulation of wild and mutant proteasome subunit beta type 8 (PSMB8) protein: Implications for restoration of inflammation in experimental autoimmune encephalomyelitis pathogenesis. *Heliyon, 11*(1), e41166. [https://doi.org/10.1016/j.heliyon.2024.e41166](https://doi.org/10.1016/j.heliyon.2024.e41166)

-----

## 📄 License

This project is licensed under the MIT License. See the `LICENSE` file for details.
