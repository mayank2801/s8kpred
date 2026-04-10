#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 16:31:05 2024

@author: mayank
"""
# The scikit-learn version is 1.5.2.
# The Pandas version is 2.1.4.
# The csv version is 1.0.
# The joblib version is 1.2.0.
# The numpy version is 1.26.4.
import sklearn
import sys
print('The scikit-learn version is {}.'.format(sklearn.__version__))

from Bio import SeqIO
import pandas as pd
print('The Pandas version is {}.'.format(pd.__version__))
import csv
print('The csv version is {}.'.format(csv.__version__))
# import joblib
# print('The joblib version is {}.'.format(joblib.__version__))
# from sklearn.preprocessing import LabelEncoder 
import numpy as np 
print('The numpy version is {}.'.format(np.__version__))
import xgboost as xgb
print('The xgboost version is {}.'.format(xgb.__version__))
# Define the path to your FASTA file

import pandas as pd
def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))

def S8kPred_format_ss2(data, ss, ss_conf,ID):
    ''' Formats output for S8kPred .ss2 files. 
    ''' 
    lines = ['# S8kPred VFORMAT \n'+str(ID)+'\n # AA  SS  B  E  G  H  I  L  S  T \n']
    for i in range(len(ss)):
        lines.append("%4d %c %c  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (i + 1, data[i], ss[i], ss_conf[i,0], ss_conf[i,1], ss_conf[i,2] , ss_conf[i,3], ss_conf[i,4], ss_conf[i,5], ss_conf[i,6], ss_conf[i,7]))
    return lines

def S8kPred_format_fas(data, ss, ss_conf, ID, include_conf=False):
    ''' Formats output as a pseudo-FASTA file
    ''' 
    lines=['> S8kPred Fasta '+str(ID) +'\n']
    lines.append("".join([d for d in data]))
    lines.append("".join([s.item() for s in ss]))
    
    if include_conf:
        lines.append(np.array2string(ss_conf[:,0],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,1],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,2],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,3],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,4],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,5],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,6],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
        lines.append(np.array2string(ss_conf[:,7],max_line_width=1e6, precision=3,formatter={'float_kind':lambda x: "%.3f" % x}).replace('[','').replace(']',''))
    return lines
    
def S8kPred_format_horiz(data, ss, ss_conf,ID):
    ''' Formats output for the S8kPred HFORMAT .horiz files. 
        Care must be taken as there is a fixed column width of 60 char
    '''    
 
    lines=['# S8kPred HFORMAT  '+str(ID) +'\n']
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
def EightStateSSPred(JobID):
    # Load the data into a pandas DataFrame
    
    
    Propensitydf=pd.read_csv('Script/Data/TriPeptidePropensityEightStateSecStructure.csv')
    Propensitydf.columns
    AminoAcidBinaryTableDf=pd.read_csv('Script/Data/TripeptideBinaryTable_60.csv')
    
    
    # fasta_file = "example.fasta"  # Replace with your file path
    fasta_file= "Jobs/"+JobID+"/FASTA/input_sequence.fasta" 
    i=1
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        Inputdf = pd.DataFrame(columns=['ID','sequence','Residue9th','PR1_B', 'PR1_E', 'PR1_G', 'PR1_H', 'PR1_I', 'PR1_L', 'PR1_S', 'PR1_T', 'PR2_B', 'PR2_E', 'PR2_G', 'PR2_H', 'PR2_I', 'PR2_L', 'PR2_S', 'PR2_T', 'PR3_B', 'PR3_E', 'PR3_G', 'PR3_H', 'PR3_I', 'PR3_L', 'PR3_S', 'PR3_T', 'PR4_B', 'PR4_E', 'PR4_G', 'PR4_H', 'PR4_I', 'PR4_L', 'PR4_S', 'PR4_T', 'PR5_B', 'PR5_E', 'PR5_G', 'PR5_H', 'PR5_I', 'PR5_L', 'PR5_S', 'PR5_T', 'PR6_B', 'PR6_E', 'PR6_G', 'PR6_H', 'PR6_I', 'PR6_L', 'PR6_S', 'PR6_T', 'PR7_B', 'PR7_E', 'PR7_G', 'PR7_H', 'PR7_I', 'PR7_L', 'PR7_S', 'PR7_T', 'PR8_B', 'PR8_E', 'PR8_G', 'PR8_H', 'PR8_I', 'PR8_L', 'PR8_S', 'PR8_T', 'PR9_B', 'PR9_E', 'PR9_G', 'PR9_H', 'PR9_I', 'PR9_L', 'PR9_S', 'PR9_T', 'PR10_B', 'PR10_E', 'PR10_G', 'PR10_H', 'PR10_I', 'PR10_L', 'PR10_S', 'PR10_T', 'PR11_B', 'PR11_E', 'PR11_G', 'PR11_H', 'PR11_I', 'PR11_L', 'PR11_S', 'PR11_T', 'PR12_B', 'PR12_E', 'PR12_G', 'PR12_H', 'PR12_I', 'PR12_L', 'PR12_S', 'PR12_T', 'PR13_B', 'PR13_E', 'PR13_G', 'PR13_H', 'PR13_I', 'PR13_L', 'PR13_S', 'PR13_T', 'PR14_B', 'PR14_E', 'PR14_G', 'PR14_H', 'PR14_I', 'PR14_L', 'PR14_S', 'PR14_T', 'PR15_B', 'PR15_E', 'PR15_G', 'PR15_H', 'PR15_I', 'PR15_L', 'PR15_S', 'PR15_T', 'PR16_B', 'PR16_E', 'PR16_G', 'PR16_H', 'PR16_I', 'PR16_L', 'PR16_S', 'PR16_T', 'PR17_B', 'PR17_E', 'PR17_G', 'PR17_H', 'PR17_I', 'PR17_L', 'PR17_S', 'PR17_T',
        'R1_A','R1_C','R1_D','R1_E','R1_F','R1_G','R1_H','R1_I','R1_K','R1_L','R1_M','R1_N','R1_P','R1_Q','R1_R','R1_S','R1_T','R1_V','R1_W','R1_Y','R2_A','R2_C','R2_D','R2_E','R2_F','R2_G','R2_H','R2_I','R2_K','R2_L','R2_M','R2_N','R2_P','R2_Q','R2_R','R2_S','R2_T','R2_V','R2_W','R2_Y','R3_A','R3_C','R3_D','R3_E','R3_F','R3_G','R3_H','R3_I','R3_K','R3_L','R3_M','R3_N','R3_P','R3_Q','R3_R','R3_S','R3_T','R3_V','R3_W','R3_Y'
        ])
        
        print(f"ID: {record.id}")
        print(f"Description: {record.description}")
        print(f"Sequence: {record.seq}\n")
        Sequence=record.seq
        ID=record.id
        SeqID='Seq_'+str(i)
        i=i+1
        print(f"Sequence Length: {len(record.seq)}\n")
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
                    Propensity.extend(Propensitydf.iloc[IndexListofPropensity[0],2:10].to_list())
                      
                    # IndexListofAminoAcidBinary = AminoAcidBinaryTableDf.index[AminoAcidBinaryTableDf['Residue'] == Sequence[k]].tolist()
                    # SequenceBinary.extend(AminoAcidBinaryTableDf.iloc[IndexListofAminoAcidBinary[0],1:21].to_list())
                
                MTriPeptide=Sequence[j+7]+Sequence[j+8]+Sequence[j+9]
                updated_MTriPeptide = MTriPeptide.replace('X', 'G')
                IndexListofAminoAcidBinary = AminoAcidBinaryTableDf.index[AminoAcidBinaryTableDf['Tripeptide'] == updated_MTriPeptide].tolist()
                SequenceBinary.extend(AminoAcidBinaryTableDf.iloc[IndexListofAminoAcidBinary[0],1:61].to_list())
                # Propensity.insert(0, df['PDBID'][i])
                Propensity.insert(0, SeqID)
                Propensity.insert(1, Sequence17W)
                # Propensity.insert(2, sec_structure[j:j+17])
                Propensity.insert(2, Residue9th)
                # Propensity.insert(4, Residue9thSecStruc)
                Propensity.extend(SequenceBinary)
                # print(len(Propensity))
                # print(Inputdf.columns)
                Inputdf.loc[ len(Inputdf) ] =Propensity 
                # writer.writerow(Propensity)
    # Inputdf.to_csv('TriPeptideAnalysis/Test/7W2H.csv')
    # # f.close()           
        PSSMDf=pd.read_csv("Jobs/"+JobID+"/PSSM_Features_ML_17W.csv")
        Inputdf_PSSM_merged_df = pd.merge(Inputdf, PSSMDf, on=["ID", "sequence"], how="inner")
        if len(Inputdf_PSSM_merged_df) >1:
            # Load the model (when needed)
            model_path = 'Script/Data/TriPeptideMachineLeraning17ResidueWindowEightPropensityBinary602NDPSSM_best_xgb_model_cpu.ubj'
            trained_model = xgb.XGBClassifier()
            trained_model.load_model(model_path)
            Input=Inputdf_PSSM_merged_df.iloc[:,3:]
            Input.columns
            y_pred = trained_model.predict(Input)
            probabilities = trained_model.predict_proba(Input)
            merged_secqunece = ''.join(Inputdf_PSSM_merged_df['Residue9th'].astype(str))
        
            # Define the mapping
            replace_map = {0: 'B', 1: 'E', 2: 'G',3:'H',4: 'I',5: 'L',6: 'S',7: 'T'}
            
            # Replace values using np.vectorize
            vectorized_replace = np.vectorize(lambda x: replace_map.get(x, x))  # Default to original if not found
            y_pred = vectorized_replace(y_pred)
            Predicted_SecStructure=''.join(y_pred)
            # print(merged_secqunece)
            # print(Predicted_SecStructure)
            # print(probabilities)
            # Create the DataFrame by combining all columns
            df = pd.DataFrame({
                'Residue': Inputdf_PSSM_merged_df['Residue9th'].to_list(),
                'Secondary Structure': y_pred
            })
            
            # Concatenate with the transposed columns
            combined_df = pd.DataFrame(probabilities, columns=['B', 'E', 'G', 'H', 'I', 'L', 'S','T'])
            # df = pd.concat([df, probabilities_data], axis=1)
            df = pd.concat([df, combined_df], axis=1)
            print(merged_secqunece)
            print(Predicted_SecStructure)
            # print(probabilities.dtype)
            S8kPred_format_ss2_lines=S8kPred_format_ss2(Inputdf_PSSM_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            
            S8kPred_format_fas_lines=S8kPred_format_fas(Inputdf_PSSM_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            S8kPred_format_horiz_lines=S8kPred_format_horiz(Inputdf_PSSM_merged_df['Residue9th'].to_list(), y_pred, probabilities,ID)
            # print(lines)
            f=open("Jobs/"+JobID+"/ResultEightState.ss2", 'a')
            for line in S8kPred_format_ss2_lines:
                f.write(line+'\n')
            f.close()
            f2=open("Jobs/"+JobID+"/ResultEightState.horiz", 'a')
            for line in S8kPred_format_horiz_lines:
                f2.write(line+'\n')
            f2.close()
            f3=open("Jobs/"+JobID+"/ResultEightState.fas", 'a')
            for line in S8kPred_format_fas_lines:
                f3.write(line+'\n')
            f3.close()
            df.to_csv("Jobs/"+JobID+"/ResultEightState.csv",mode='a')
            del trained_model
            # Display the modified array
    return 1
# JobID='test'
# EightStateSSPred(JobID)
if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2])
