# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 18:18:31 2024

@author: mayank
"""

# import os
# from Bio.PDB.MMCIFParser import MMCIFParser
# from Bio.PDB.PDBIO import PDBIO

# def convert_all_mmcif_to_pdb(JobID):

#     input_folder = "Jobs/"+JobID+"/mmCIF"
#     output_folder ="Jobs/"+JobID+"/PDB"
#     if not os.path.exists(output_folder):
#         os.makedirs(output_folder)  # Create the output folder if it doesn't exist

#     parser = MMCIFParser(QUIET=True)
#     io = PDBIO()

#     for filename in os.listdir(input_folder):
#         if filename.endswith(".cif"):
#             input_path = os.path.join(input_folder, filename)
#             output_path = os.path.join(output_folder, filename.replace(".cif", ".pdb"))
            
#             try:
#                 structure = parser.get_structure(filename, input_path)
#                 io.set_structure(structure)
#                 io.save(output_path)
#                 print(f"Converted: {filename} -> {output_path}")
#             except Exception as e:
#                 print(f"Error converting {filename}: {e}")

# JobID='RamPlot675bdb5ff3557'
# convert_all_mmcif_to_pdb(JobID)
import os
import sys
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from datetime import datetime

def convert_all_mmcif_to_pdb(JobID):
    """
    Converts all MMCIF files in a JobID folder to PDB format.
    Logs the conversion status and returns 1 for success or 0 for failure.

    Parameters:
        JobID (str): Job identifier for the input and output folders.

    Returns:
        int: 1 if successful, 0 otherwise.
    """
    input_folder = f"Jobs/{JobID}/mmCIF"
    output_folder = f"Jobs/{JobID}/PDB"
    log_file = f"Jobs/{JobID}/log.dat"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # Create the output folder if it doesn't exist

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    success = True  # Flag to track overall success

    # Open the log file
    with open(log_file, "a") as log:
        log.write(f"\n--- Conversion Log: {datetime.now()} ---\n")

        for filename in os.listdir(input_folder):
            if filename.endswith(".cif"):
                input_path = os.path.join(input_folder, filename)
                output_path = os.path.join(output_folder, filename.replace(".cif", ".pdb"))

                try:
                    structure = parser.get_structure(filename, input_path)
                    io.set_structure(structure)
                    io.save(output_path)
                    log.write(f"SUCCESS: {filename} converted to {output_path}\n")
                    print(f"Converted: {filename} -> {output_path}")
                except Exception as e:
                    log.write(f"ERROR: Failed to convert {filename} - {e}\n")
                    print(f"Error converting {filename}: {e}")
                    success = False  # Set flag to False if any file fails

    # Return 1 for success, 0 for failure
    log.close()
    return 1 if success else 0


# Example usage
# JobID = 'RamPlot675bdb5ff3557'
# status = convert_all_mmcif_to_pdb(JobID)
# print(f"Conversion status: {'Success' if status == 1 else 'Failure'}")
if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])
