import os
import sys
import csv
from Bio import PDB

def run_dssp(JobID):
    """
    Run DSSP on a given PDB or CIF file and output results into a CSV file.

    Parameters:
    input_file (str): Path to PDB or mmCIF file.
    output_csv (str): Path to save the output CSV file.
    """
    
    files = os.listdir('Jobs/'+JobID+'/PDB/')
    dssp_data = []
    # dssp_output_all=''
    for pdb_code in files:      
        #log.write(pdb_code+"\n")
        print(pdb_code)
        PDBID=pdb_code.split(".")[0]
        #PDBID='ABC'
        input_file = 'Jobs/'+JobID+'/PDB/'+pdb_code
        # Load structure
        parser = PDB.PDBParser(QUIET=True) if input_file.endswith(".pdb") else PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure("protein", input_file)

        # Use the first model
        model = structure[0]

        # Run DSSP
        dssp = PDB.DSSP(model, input_file, acc_array='Sander')
        # dssp_output_all=dssp_output_all+str(dssp)
        # Extract DSSP keys and count
        dssp_keys = list(dssp.keys())
        total_residues = len(dssp)

        # Extract DSSP data
        
        for i in range(0, total_residues):
            try:
                current_key = dssp_keys[i]

                # Extract DSSP parameters
                residues = dssp[current_key][1]  # Residue type
                residue_num = current_key[1][1]  # Residue number
                chain = current_key[0]  # Chain ID
                rsa = round(dssp[current_key][3],3) # Relative Solvent Accessibility
                phi = dssp[current_key][4]  # Phi angle
                psi = dssp[current_key][5]  # Psi angle
                
                # Additional DSSP parameters
                NH_O1 = dssp[current_key][6]  # NH->O1 hydrogen bond partner
                NH_O1_E = dssp[current_key][7]  # NH->O1 hydrogen bond energy
                O_NH1 = dssp[current_key][8]  # O->NH1 hydrogen bond partner
                O_NH1_E = dssp[current_key][9]  # O->NH1 hydrogen bond energy
                NH_O2 = dssp[current_key][10]  # NH->O2 hydrogen bond partner
                NH_O2_E = dssp[current_key][11]  # NH->O2 hydrogen bond energy
                O_NH2 = dssp[current_key][12]  # O->NH2 hydrogen bond partner
                O_NH2_E = dssp[current_key][13]  # O->NH2 hydrogen bond energy


                # Define secondary structure: Use 'L' for loop if '-'
                sec_structure = 'L' if dssp[current_key][2] == '-' else dssp[current_key][2]

                # Append extracted data
                dssp_data.append([
                    PDBID,residue_num, chain, residues, sec_structure, rsa, phi, psi,
                    NH_O1, NH_O1_E, O_NH1, O_NH1_E,
                    NH_O2, NH_O2_E, O_NH2, O_NH2_E
                ])
            except Exception as e:
                print(f"Skipping residue {i} due to error: {e}")

    # Write to CSV
    with open('Jobs/'+JobID+'/Dssp_Output.csv', "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDBID",
            "Residue_No", "Chain", "Amino Acid", "Secondary Structure", "Solvent Accessibility",
            "Phi Angle", "Psi Angle",
            "NH->O1 Partner", "NH->O1 Energy", "O->NH1 Partner", "O->NH1 Energy",
            "NH->O2 Partner", "NH->O2 Energy", "O->NH2 Partner", "O->NH2 Energy"
        ])
        writer.writerows(dssp_data)
    # with open('Jobs/'+JobID+'/Dssp_Output.txt', 'w') as txt_file:
    #     txt_file.write(dssp_output_all)

    print(f"Results saved to {JobID} Dssp_Output.csv")


if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])