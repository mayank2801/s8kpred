#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 18:05:27 2024

@author: mayank
"""
# The scikit-learn version is 1.5.2.
# The Pandas version is 2.1.4.
# The csv version is 1.0.
# The joblib version is 1.2.0.
# The numpy version is 1.26.4.
import sys
import numpy as np 
print('The numpy version is {}.'.format(np.__version__))
# =============================================================================
# Output helpers
# =============================================================================

ind2char={0:"C", 1:"H", 2:"E"}

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def S8kPred_format_ss2(data, ss, ss_conf,ID):
    ''' Formats output for S8kPred .ss2 files. 
    ''' 
    lines = ['# S8kPred VFORMAT '+str(ID) +'\n']
    for i in range(len(ss)):
        lines.append("%4d %c %c  %6.3f %6.3f %6.3f" % (i + 1, data[i], ss[i], ss_conf[i,0], ss_conf[i,1], ss_conf[i,2]))
    return lines

def S8kPred_format_fas(data, ss, ss_conf, ID,include_conf=False):
    ''' Formats output as a pseudo-FASTA file
    ''' 
    lines=['> S8kPred Fasta '+str(ID) +'\n']
    lines.append("".join([d for d in data]))
    lines.append("".join([s.item() for s in ss]))
    
    if include_conf:
        lines.append(np.array2string(ss_conf[:,0],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,1],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,2],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
    
    return lines
    
def S8kPred_format_horiz(data, ss, ss_conf,ID):
    ''' Formats output for the S8kPred HFORMAT .horiz files. 
        Care must be taken as there is a fixed column width of 60 char
    '''    
 
    lines=['# S8kPred HFORMAT '+str(ID) +'\n']
    sub_seqs =  list(chunkstring("".join([d for d in data]),60))
    # list(chunkstring(data,60))
    sub_ss   = list(chunkstring("".join([s.item() for s in ss]),60))
    
    num_len  =  int(np.floor(len(data)/10))
    num_seq  = ''.join(f'{str((i+1)*10):>10}' for i in range(num_len+1))
    num_seq  = list(chunkstring(num_seq,60))
        
    # get confidences then floor them and convert to string 
    conf_idxs = ss_conf.argmax(-1)
    confs = ss_conf[np.arange(len(conf_idxs)),conf_idxs[:]]
    confs = "".join([str(x) for x in np.floor(confs*10).astype(np.int32)])
    confs = list(chunkstring(confs,60))

    for idx, subsq in enumerate(sub_seqs):
        lines.append(f'\nConf: {confs[idx]}')
        lines.append(f'Pred: {sub_ss[idx]}')
        lines.append(f'  AA: {subsq}')
        lines.append(f'      {num_seq[idx]}\n')
        
    return lines

def ThreeStateSSPred(JobID):
    import sklearn
    print('The scikit-learn version is {}.'.format(sklearn.__version__))
    
    from Bio import SeqIO
    import pandas as pd
    print('The Pandas version is {}.'.format(pd.__version__))
    import csv
    print('The csv version is {}.'.format(csv.__version__))
    # import joblib
    # print('The joblib version is {}.'.format(joblib.__version__))
    # from sklearn.preprocessing import LabelEncoder 
    import xgboost as xgb
    # Define the path to your FASTA file
    Propensitydf=pd.read_csv('Script/Data/TriPeptidePropensityThreeStateSecStructure2AND.csv')
    Propensitydf.columns
    # AminoAcidBinaryTableDf=pd.read_csv('Script/Data/AminoAcidBinaryTable.csv')
    
    
    
    fasta_file = "Jobs/"+JobID+"/FASTA/input_sequence.fasta"  # Replace with your file path
    i=1
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        Inputdf = pd.DataFrame(columns=['ID','sequence','Residue9th','H1', 'E1', 'L1',
         'H2', 'E2', 'L2',
         'H3', 'E3', 'L3',
         'H4', 'E4', 'L4',
         'H5', 'E5', 'L5',
         'H6', 'E6', 'L6',
         'H7', 'E7', 'L7',
         'H8', 'E8', 'L8',
         'H9', 'E9', 'L9',
         'H10', 'E10', 'L10',
         'H11', 'E11', 'L11',
         'H12', 'E12', 'L12',
         'H13', 'E13', 'L13',
         'H14', 'E14', 'L14',
         'H15', 'E15', 'L15',
         'H16', 'E16', 'L16',
         'H17', 'E17', 'L17'])
        print(f"ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Sequence: {record.seq}\n")
        Sequence=record.seq
        print(f"Sequence Length: {len(record.seq)}\n")
        ID=record.id
        SeqID='Seq_'+str(i)
        i=i+1
        PrefixSuffixSequence='GGGGGGGG'
        Sequence=PrefixSuffixSequence+Sequence+PrefixSuffixSequence
        Sequence=Sequence.upper()
        SequnceLen=len(Sequence)
        for j in range(SequnceLen):
    
            Sequence17W=Sequence[j:j+17]  
            windowsize=len(Sequence17W)
            if windowsize==17:
    
                # print(Sequence17W)
                # print(sec_structure[j:j+17])
                Residue9th=Sequence[j+8]
                # print(Residue9th)
                # Residue9thSecStruc=sec_structure[j+8]
                # print(Residue9thSecStruc)
                Propensity=[]
                SequenceBinary=[]
                for k in range(j,j+17):
                    # print(k)
                    if k==0:
                        S1='G'
                        S2=Sequence[k]
                        S3=Sequence[k+1]
                    elif k+1 >= SequnceLen:
                        S1=Sequence[k-1]
                        S2=Sequence[k]
                        S3='G'
                    else:
                        S1=Sequence[k-1]
                        S2=Sequence[k]
                        S3=Sequence[k+1]
                    tripettide=S1+S2+S3
                    # print(tripettide)
                    updated_tripettide = tripettide.replace('X', 'G')
                    IndexListofPropensity = Propensitydf.index[Propensitydf['TriPeptide'] == updated_tripettide].tolist()
                    Propensity.extend(Propensitydf.iloc[IndexListofPropensity[0],2:5].to_list())
                      
                    # IndexListofAminoAcidBinary = AminoAcidBinaryTableDf.index[AminoAcidBinaryTableDf['Residue'] == Sequence[k]].tolist()
                    # SequenceBinary.extend(AminoAcidBinaryTableDf.iloc[IndexListofAminoAcidBinary[0],1:21].to_list())
                # Propensity.insert(0, df['PDBID'][i])
                Propensity.insert(0, SeqID)
                Propensity.insert(1, Sequence17W)
                # Propensity.insert(2, sec_structure[j:j+17])
                Propensity.insert(2, Residue9th)
                # Propensity.insert(4, Residue9thSecStruc)
                # Propensity.extend(SequenceBinary)
                # print(Propensity)
                Inputdf.loc[ len(Inputdf) ] =Propensity
                
                # writer.writerow(Propensity)
    # Inputdf.to_csv('TriPeptideAnalysis/Test/7W2H.csv')
    # # f.close()           
        PSSMDf=pd.read_csv("Jobs/"+JobID+"/PSSM_Features_ML_17W.csv")
        BinaryDf=pd.read_csv('Script/Data/TripeptideBinaryTable_60.csv')
        # Merge the DataFrames on the first two columns (assuming they are 'PDBID' and 'sequence')
        Inputdf_PSSM_merged_df = pd.merge(Inputdf, PSSMDf, on=["ID", "sequence"], how="inner")
        # Extract 8th, 9th, and 10th residues as tripeptides
        # Inputdf_PSSM_merged_df['Tripeptide'] = Inputdf_PSSM_merged_df['sequence'].str[7:10]  # 0-based indexing, extract characters 8-10
        # Inputdf_PSSM_merged_df['Tripeptide'] = Inputdf_PSSM_merged_df['sequence'].str.replace('X', 'G').str[7:10]
        Inputdf_PSSM_merged_df['Tripeptide'] = (
    Inputdf_PSSM_merged_df['sequence']
    .astype(str)  # Convert to string to avoid AttributeError
    .str.replace('X', 'G')  # Replace 'X' with 'G'
    .str[7:10]  # Extract tripeptide from positions 7 to 9
)

        # Merge with binary table to get binary values for each tripeptide
        Inputdf_PSSM_Binary_merged_df = pd.merge(Inputdf_PSSM_merged_df, BinaryDf, how='left', left_on='Tripeptide', right_on='Tripeptide')
        # updated_df = pd.merge(df, df2, how='left', left_on='Residue9th', right_on='AminoAcid')
    
        
        # Save the updated DataFrame to a new CSV file
        Inputdf_PSSM_Binary_merged_df  = Inputdf_PSSM_Binary_merged_df.drop(columns=['Tripeptide'])
        # Remove rows with NaN values (optional, based on your requirements)
        Inputdf_PSSM_Binary_merged_df  = Inputdf_PSSM_Binary_merged_df.dropna()
        if len(Inputdf_PSSM_Binary_merged_df) >1:
            # Load the model (when needed)
            model_path = 'Script/Data/TriPeptideML17WPropensity2ANDThreeStateWithBinary60PssmMerged_best_xgb_model_cpu.json'
            trained_model = xgb.XGBClassifier()
            trained_model.load_model(model_path)
            Input=Inputdf_PSSM_Binary_merged_df.iloc[:,3:]
            Input.columns
            y_pred = trained_model.predict(Input)
            probabilities = trained_model.predict_proba(Input)
            merged_secqunece = ''.join(Inputdf_PSSM_Binary_merged_df['Residue9th'].astype(str))
        
            # Define the mapping
            replace_map = {0: 'E', 1: 'H', 2: 'L'}
            
            # Replace values using np.vectorize
            vectorized_replace = np.vectorize(lambda x: replace_map.get(x, x))  # Default to original if not found
            y_pred = vectorized_replace(y_pred)
            Predicted_SecStructure=''.join(y_pred)
            # Inputdf_PSSM_Binary_merged_df['SecStructure']=y_pred
            # Transpose the combined data to match the shape for columns
            # probabilities_data = pd.DataFrame(zip(*probabilities), columns=['Col3', 'Col4', 'Col5'])
            
            # Create the DataFrame by combining all columns
            df = pd.DataFrame({
                'Residue': Inputdf_PSSM_Binary_merged_df['Residue9th'].to_list(),
                'Secondary Structure': y_pred,'ID': Inputdf_PSSM_Binary_merged_df['ID'].to_list(),
            })
            
            # Concatenate with the transposed columns
            combined_df = pd.DataFrame(probabilities, columns=['E', 'H', 'L'])
            # df = pd.concat([df, probabilities_data], axis=1)
            df = pd.concat([df, combined_df], axis=1)
            print(merged_secqunece)
            print(Predicted_SecStructure)
            # print(probabilities.dtype)
            S8kPred_format_ss2_lines=S8kPred_format_ss2(Inputdf_PSSM_Binary_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            
            S8kPred_format_fas_lines=S8kPred_format_fas(Inputdf_PSSM_Binary_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            S8kPred_format_horiz_lines=S8kPred_format_horiz(Inputdf_PSSM_Binary_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            # print(lines)
            f=open("Jobs/"+JobID+"/ResultThreeState.ss2", 'a')
            for line in S8kPred_format_ss2_lines:
                f.write(line+'\n')
            f.close()
            f2=open("Jobs/"+JobID+"/ResultThreeState.horiz", 'a')
            for line in S8kPred_format_horiz_lines:
                f2.write(line+'\n')
            f2.close()
            f3=open("Jobs/"+JobID+"/ResultThreeState.fas", 'a')
            for line in S8kPred_format_fas_lines:
                f3.write(line+'\n')
            f3.close()
            df.to_csv("Jobs/"+JobID+"/ResultThreeState.csv",mode='a')
            del trained_model
            # Generate Catoon Plot
            import subprocess

            script_path = "Script/PlotSecondaryStructureCartoon.py"
            # conda_python = "C:/Users/mayank/AppData/Local/Schrodinger/PyMOL2/envs/biotite/python"
            conda_python = "/var/www/miniconda3/envs/biotite/bin/python"
            
            fun = 'PlotSecondaryStructureCartoon'
            result = subprocess.run([conda_python, script_path, fun, Predicted_SecStructure, JobID,SeqID], capture_output=True, text=True)
        # Display the modified array
    return 1
# JobID='test'
# ThreeStateSSPred(JobID)
if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])
    