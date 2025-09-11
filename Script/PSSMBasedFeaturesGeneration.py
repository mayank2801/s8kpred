#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 11:07:33 2025

@author: mayank
"""

from Bio import SeqIO
import subprocess
import os
import time
import multiprocessing
import re
import sys
import numpy as np
import csv

def read_pssm_matrix(input_matrix):
    """Reads a PSSM file and extracts the amino acid sequence and numerical matrix."""
    PSSM = []
    amino_acids = []  # List to store the amino acid sequence
    numeric_pattern = re.compile(r'-?\d+')

    with open(input_matrix, 'r') as stream:
        header_skipped = False  # To ensure the header is skipped
        for line in stream:
            if not header_skipped:
                # Skip until we find the matrix start (the first line after the header)
                if line.startswith("Last position-specific scoring matrix computed"):
                    header_skipped = True
                continue
            
            columns = line.split()
            if not columns or not columns[0].isdigit():
                # Stop processing when reaching non-matrix content
                continue
            
            # Extract amino acid only from valid PSSM lines
            if len(columns) >= 2 and columns[0].isdigit():
                amino_acids.append(columns[1])  # Second column contains the amino acid
            
            # Extract numeric values for the PSSM matrix
            if len(columns) >= 22:  # Ensure enough columns for numerical data
                values = [int(x) for x in columns[2:22]]
                PSSM.append(values)

    return ''.join(amino_acids), np.array(PSSM)

def ensure_csv_exists(output_file, header):
    """
    Ensure that the CSV file exists. If not, create it with the provided header.
    """
    file_exists = os.path.isfile(output_file)
    
    with open(output_file, mode="a", newline="") as out_file:
        csv_writer = csv.writer(out_file)
        
        # Write the header only if the file is newly created
        if not file_exists:
            csv_writer.writerow(header)
def generate_pssm_files(JobID):
    # Define input and output parameters
    # JobID='test'
    input_fasta="Jobs/"+JobID+"/FASTA/input_sequence.fasta" 
    # input_fasta = "test2.fasta"
    database = "UNI_REF_Protein_Database_NCBI/uniref50"
    num_iterations = 3
    output_dir = "Jobs/"+JobID+"/pssm_outputs"
    num_threads=4
    # Get the number of available CPU cores for multi-threading
    if num_threads=='':
        num_threads = multiprocessing.cpu_count()
        
    
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Open a log file to record time taken for each sequence
    log_file = os.path.join("Jobs/"+JobID, "log.dat")
    with open(log_file, "a") as log:
        log.write("PSSM Calculation\n")
        log.write("Sequence_ID\tTime_Taken(s)\tStatus\n")  # Header for log file
        i=1
        # Process each sequence in the multi-FASTA file
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_id = record.id  # Get sequence ID
            SeqID='Seq_'+str(i)
            i=i+1
            temp_fasta = f"{output_dir}/{SeqID}.fasta"  # Temporary FASTA file
            pssm_output = os.path.join(output_dir, f"{SeqID}.pssm")  # Output PSSM file
    
            # Write individual sequence to a temporary FASTA file
            with open(temp_fasta, "w") as f:
                f.write(f">{record.id}\n{record.seq}\n")
    
            print(f"Processing {seq_id} with {num_threads} threads...")
    
            # Measure execution time
            start_time = time.time()
            
            # Run PSI-BLAST command
            cmd = [
                "Script/ncbi-blast-2.13.0/bin/psiblast",
                "-query", temp_fasta,
                "-db", database,
                "-num_iterations", str(num_iterations),
                "-num_threads", str(num_threads),
                "-out_ascii_pssm", pssm_output
            ]
    
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                time_taken = round(time.time() - start_time, 2)  # Calculate time elapsed
    
                # Check if PSSM file was generated
                if not os.path.exists(pssm_output) or os.path.getsize(pssm_output) == 0:
                    raise FileNotFoundError(f"PSSM file missing for {seq_id}")
    
                print(f"Completed {seq_id} in {time_taken} seconds.")
                log.write(f"{seq_id}\t{time_taken}\tSUCCESS\n")
    
            except subprocess.CalledProcessError:
                print(f"Error: PSI-BLAST failed for {seq_id}")
                log.write(f"{seq_id}\tERROR\tPSI-BLAST Failed\n")
            except FileNotFoundError as e:
                print(f"Error: {e}")
                log.write(f"{seq_id}\tERROR\tPSSM Not Generated\n")
            except Exception as e:
                print(f"Unexpected error for {seq_id}: {e}")
                log.write(f"{seq_id}\tERROR\t{e}\n")
    
            # Remove temporary FASTA file
            # if temp_fasta!=input_fasta:
            os.remove(temp_fasta)
    
    print("PSSM generation completed! Check the 'pssm_outputs' folder and log file.")







# def process_pssm_files(JobID):
    """Processes all PSSM files in a directory and writes scaled features with sequences to a CSV file."""
    # files = os.listdir(input_dir)
    input_dir = "Jobs/"+JobID+"/pssm_outputs"
    files = [f for f in os.listdir(input_dir) if f.endswith(".pssm")]
    if files!='':
        log=open("Jobs/"+JobID+"/log.dat",mode='a')
        output_file="Jobs/"+JobID+"/PSSM_Features_ML_17W.csv"
        header=['ID','sequence','R1_P_A', 'R1_P_R', 'R1_P_N', 'R1_P_D', 'R1_P_C', 'R1_P_Q', 'R1_P_E', 'R1_P_G', 'R1_P_H', 'R1_P_I', 'R1_P_L', 'R1_P_K', 'R1_P_M', 'R1_P_F', 'R1_P_P', 'R1_P_S', 'R1_P_T', 'R1_P_W', 'R1_P_Y', 'R1_P_V', 'R2_P_A', 'R2_P_R', 'R2_P_N', 'R2_P_D', 'R2_P_C', 'R2_P_Q', 'R2_P_E', 'R2_P_G', 'R2_P_H', 'R2_P_I', 'R2_P_L', 'R2_P_K', 'R2_P_M', 'R2_P_F', 'R2_P_P', 'R2_P_S', 'R2_P_T', 'R2_P_W', 'R2_P_Y', 'R2_P_V', 'R3_P_A', 'R3_P_R', 'R3_P_N', 'R3_P_D', 'R3_P_C', 'R3_P_Q', 'R3_P_E', 'R3_P_G', 'R3_P_H', 'R3_P_I', 'R3_P_L', 'R3_P_K', 'R3_P_M', 'R3_P_F', 'R3_P_P', 'R3_P_S', 'R3_P_T', 'R3_P_W', 'R3_P_Y', 'R3_P_V', 'R4_P_A', 'R4_P_R', 'R4_P_N', 'R4_P_D', 'R4_P_C', 'R4_P_Q', 'R4_P_E', 'R4_P_G', 'R4_P_H', 'R4_P_I', 'R4_P_L', 'R4_P_K', 'R4_P_M', 'R4_P_F', 'R4_P_P', 'R4_P_S', 'R4_P_T', 'R4_P_W', 'R4_P_Y', 'R4_P_V', 'R5_P_A', 'R5_P_R', 'R5_P_N', 'R5_P_D', 'R5_P_C', 'R5_P_Q', 'R5_P_E', 'R5_P_G', 'R5_P_H', 'R5_P_I', 'R5_P_L', 'R5_P_K', 'R5_P_M', 'R5_P_F', 'R5_P_P', 'R5_P_S', 'R5_P_T', 'R5_P_W', 'R5_P_Y', 'R5_P_V', 'R6_P_A', 'R6_P_R', 'R6_P_N', 'R6_P_D', 'R6_P_C', 'R6_P_Q', 'R6_P_E', 'R6_P_G', 'R6_P_H', 'R6_P_I', 'R6_P_L', 'R6_P_K', 'R6_P_M', 'R6_P_F', 'R6_P_P', 'R6_P_S', 'R6_P_T', 'R6_P_W', 'R6_P_Y', 'R6_P_V', 'R7_P_A', 'R7_P_R', 'R7_P_N', 'R7_P_D', 'R7_P_C', 'R7_P_Q', 'R7_P_E', 'R7_P_G', 'R7_P_H', 'R7_P_I', 'R7_P_L', 'R7_P_K', 'R7_P_M', 'R7_P_F', 'R7_P_P', 'R7_P_S', 'R7_P_T', 'R7_P_W', 'R7_P_Y', 'R7_P_V', 'R8_P_A', 'R8_P_R', 'R8_P_N', 'R8_P_D', 'R8_P_C', 'R8_P_Q', 'R8_P_E', 'R8_P_G', 'R8_P_H', 'R8_P_I', 'R8_P_L', 'R8_P_K', 'R8_P_M', 'R8_P_F', 'R8_P_P', 'R8_P_S', 'R8_P_T', 'R8_P_W', 'R8_P_Y', 'R8_P_V', 'R9_P_A', 'R9_P_R', 'R9_P_N', 'R9_P_D', 'R9_P_C', 'R9_P_Q', 'R9_P_E', 'R9_P_G', 'R9_P_H', 'R9_P_I', 'R9_P_L', 'R9_P_K', 'R9_P_M', 'R9_P_F', 'R9_P_P', 'R9_P_S', 'R9_P_T', 'R9_P_W', 'R9_P_Y', 'R9_P_V', 'R10_P_A', 'R10_P_R', 'R10_P_N', 'R10_P_D', 'R10_P_C', 'R10_P_Q', 'R10_P_E', 'R10_P_G', 'R10_P_H', 'R10_P_I', 'R10_P_L', 'R10_P_K', 'R10_P_M', 'R10_P_F', 'R10_P_P', 'R10_P_S', 'R10_P_T', 'R10_P_W', 'R10_P_Y', 'R10_P_V', 'R11_P_A', 'R11_P_R', 'R11_P_N', 'R11_P_D', 'R11_P_C', 'R11_P_Q', 'R11_P_E', 'R11_P_G', 'R11_P_H', 'R11_P_I', 'R11_P_L', 'R11_P_K', 'R11_P_M', 'R11_P_F', 'R11_P_P', 'R11_P_S', 'R11_P_T', 'R11_P_W', 'R11_P_Y', 'R11_P_V', 'R12_P_A', 'R12_P_R', 'R12_P_N', 'R12_P_D', 'R12_P_C', 'R12_P_Q', 'R12_P_E', 'R12_P_G', 'R12_P_H', 'R12_P_I', 'R12_P_L', 'R12_P_K', 'R12_P_M', 'R12_P_F', 'R12_P_P', 'R12_P_S', 'R12_P_T', 'R12_P_W', 'R12_P_Y', 'R12_P_V', 'R13_P_A', 'R13_P_R', 'R13_P_N', 'R13_P_D', 'R13_P_C', 'R13_P_Q', 'R13_P_E', 'R13_P_G', 'R13_P_H', 'R13_P_I', 'R13_P_L', 'R13_P_K', 'R13_P_M', 'R13_P_F', 'R13_P_P', 'R13_P_S', 'R13_P_T', 'R13_P_W', 'R13_P_Y', 'R13_P_V', 'R14_P_A', 'R14_P_R', 'R14_P_N', 'R14_P_D', 'R14_P_C', 'R14_P_Q', 'R14_P_E', 'R14_P_G', 'R14_P_H', 'R14_P_I', 'R14_P_L', 'R14_P_K', 'R14_P_M', 'R14_P_F', 'R14_P_P', 'R14_P_S', 'R14_P_T', 'R14_P_W', 'R14_P_Y', 'R14_P_V', 'R15_P_A', 'R15_P_R', 'R15_P_N', 'R15_P_D', 'R15_P_C', 'R15_P_Q', 'R15_P_E', 'R15_P_G', 'R15_P_H', 'R15_P_I', 'R15_P_L', 'R15_P_K', 'R15_P_M', 'R15_P_F', 'R15_P_P', 'R15_P_S', 'R15_P_T', 'R15_P_W', 'R15_P_Y', 'R15_P_V', 'R16_P_A', 'R16_P_R', 'R16_P_N', 'R16_P_D', 'R16_P_C', 'R16_P_Q', 'R16_P_E', 'R16_P_G', 'R16_P_H', 'R16_P_I', 'R16_P_L', 'R16_P_K', 'R16_P_M', 'R16_P_F', 'R16_P_P', 'R16_P_S', 'R16_P_T', 'R16_P_W', 'R16_P_Y', 'R16_P_V', 'R17_P_A', 'R17_P_R', 'R17_P_N', 'R17_P_D', 'R17_P_C', 'R17_P_Q', 'R17_P_E', 'R17_P_G', 'R17_P_H', 'R17_P_I', 'R17_P_L', 'R17_P_K', 'R17_P_M', 'R17_P_F', 'R17_P_P', 'R17_P_S', 'R17_P_T', 'R17_P_W', 'R17_P_Y', 'R17_P_V']
        # Ensure the CSV file exists with the header
        ensure_csv_exists(output_file, header)
        with open(output_file, mode="a", newline="") as out_file:
            csv_writer = csv.writer(out_file)
            # csv_writer.writerow(header)
            for file in files:
                try:
                    file_path = os.path.join(input_dir, file)
                    ID = os.path.splitext(file)[0]
                    print(f"Processing {file}...")
                    
                    # Extract the amino acid sequence and PSSM matrix
                    amino_acid_sequence, PSSM = read_pssm_matrix(file_path)
                    selected_columns = PSSM[:, :20].astype(float)  # Select first 20 columns
                    scaled_array = 1 / (1 + np.exp(-selected_columns))  # Apply logistic function
                    flattened_array =  scaled_array.flatten(order='C')
                    
                    # Insert ID and amino acid sequence at the start of the row
                    # flattened_array = flattened_array.tolist()
                    flattened_array = [round(x, 4) for x in flattened_array.tolist()]
                    PrefixSuffixSequence='GGGGGGGG'
                    Sequence=PrefixSuffixSequence+amino_acid_sequence+PrefixSuffixSequence
                    # Generate 60 random values in the range [0, 0.3]
                    random_values = np.round(np.random.uniform(0, 0.3, 160),4) # Shape (60, 20)
                    
                    # Concatenate at the start and end
                    flattened_array=np.concatenate((random_values, flattened_array, random_values))
                    # new_array = np.vstack((random_values, flattened_array, random_values))
                    SequnceLen=len(Sequence)
                    # Sequence=amino_acid_sequence
                    
                    for j in range(SequnceLen):
                        OutPutRow=[]
                        Sequence17W=Sequence[j:j+17]  
                        windowsize=len(Sequence17W)
                        if windowsize==17:
                            
                            # print(Sequence17W)
                            # print(sec_structure[j:j+17])
                            # Residue9th=Sequence[j+9]
                            start=j*20
                            OutPutRow.extend(flattened_array[start:start+340])
                            OutPutRow.insert(0, ID)
                            OutPutRow.insert(1, Sequence17W)
                            csv_writer.writerow(OutPutRow)
                    
                except Exception as e:
                    print(f"Error processing {file}: {e}")
                    log.write(f"Error processing {file}: {e}")
    else:
        print('No PSSM File Generated')
        log.write("No PSSM File Generated")
    log.close()
# Paths
# input_directory = "Jobs/"+JobID+"/pssm_outputs"
# output_csv = "Jobs/"+JobID+"/PSSM_Features_ML_17W.csv"
# Process files
# JobID='test'
# generate_pssm_files(JobID)
if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])