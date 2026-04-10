# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 14:54:50 2025

@author: mayank
"""


import os
import sys
import csv
import subprocess

def parse_stride_output(stride_output, pdb_id):
    """Parse STRIDE output into a structured format."""
    data = []
    for line in stride_output.splitlines():
        if line.startswith("ASG"):
            parts = line.strip().split()
            if len(parts) >= 11:
                data.append({
                    'PDB_ID': pdb_id,
                    'ResNum': parts[3],
                    'Residue': parts[1],
                    'Chain': parts[2],
                    'Code': parts[5],
                    'Structure': parts[6],
                    'Solvent Accessibility': float(parts[9]),
                    'Phi': float(parts[7]),
                    'Psi': float(parts[8])
                    
                    
                })
    return data

def write_outputs(data, stride_output_all,txt_path, csv_path):
    """Write parsed data to text and CSV files."""
    with open(txt_path, 'w') as txt_file:
        txt_file.write(stride_output_all)
        # for row in data:
        #     txt_file.write(f"{row['Residue']:>4} {row['ResNum']:>4} {row['Structure']:>10} "
        #                    f"Phi={row['Phi']:>7.2f}, Psi={row['Psi']:>7.2f}, "
        #                    f"Area={row['Solvent Accessibility']:>6.1f}, PDB={row['PDB_ID']}\n")

    with open(csv_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=data[0].keys())
        writer.writeheader()
        writer.writerows(data)

    print(f"\n Parsed {len(data)} total entries.\n📁 Output:\n  - {txt_path}\n  - {csv_path}")

def run_stride(JobID):
    """Main function to run STRIDE on all PDBs in the job directory."""
    pdb_dir = os.path.join("Jobs", JobID, "PDB")
    stride_output_all=''
    if not os.path.isdir(pdb_dir):
        print(f"[ERROR] Directory not found: {pdb_dir}")
        sys.exit(1)

    pdb_files = [f for f in os.listdir(pdb_dir) if f.lower().endswith('.pdb')]
    print(f" Found {len(pdb_files)} PDB file(s) in {pdb_dir}\n")

    all_data = []

    for pdb_file in pdb_files:
        full_path = os.path.join(pdb_dir, pdb_file)

        try:
            result = subprocess.run(['./Script/stride/stride', full_path],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True,
                                    check=True)
            stride_output = result.stdout
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] STRIDE failed on {full_path}:\n{e.stderr}")
            continue  # Skip to next file

        if stride_output:
            pdb_id = os.path.splitext(pdb_file)[0]
            data = parse_stride_output(stride_output, pdb_id)
            stride_output_all=stride_output_all+stride_output
            all_data.extend(data)

    if all_data:
        out_dir = os.path.join("Jobs", JobID)
        write_outputs(all_data,stride_output_all,
                      txt_path=os.path.join(out_dir, "stride_output.txt"),
                      csv_path=os.path.join(out_dir, "stride_output.csv"))
    else:
        print("[WARNING] No valid STRIDE data parsed.")

if __name__ == "__main__":
    # Usage: python script.py run_stride <JobID>
    if len(sys.argv) != 3:
        print("Usage: python script.py run_stride <JobID>")
        sys.exit(1)

    func_name = sys.argv[1]
    job_id = sys.argv[2]

    if func_name in globals():
        globals()[func_name](job_id)
    else:
        print(f"[ERROR] Function '{func_name}' not found.")
