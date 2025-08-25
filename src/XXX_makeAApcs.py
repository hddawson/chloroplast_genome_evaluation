from aaindex import aaindex1
import pandas as pd
import os
import glob
from Bio import SeqIO
import csv
from collections import defaultdict
from tqdm import tqdm
import numpy as np

def download_all_idxs_to_csv(outfile):
    codes = aaindex1.record_codes()
    indeces = {}
    for code in codes:
        indeces[code] = aaindex1[code].values
        #'values': {'-': 0, 'A': 0.7, 'C': 0.65, 'D': 0.98, 'E': 1.04, 'F': 0.93, 'G': 1.41, 'H': 1.22, 'I': 0.78, 'K': 1.01, 'L': 0.85, 'M': 0.83, 'N': 1.42, 'P': 1.1, 'Q': 0.75, 'R': 0.34, 'S': 1.55, 'T': 1.09, 'V': 0.75, 'W': 0.62, 'Y': 0.99}

    print(indeces)
    df = pd.DataFrame(indeces)
    #remove the "-" row
    df = df.drop("-", axis=0)
    df.to_csv(outfile)

def download_idx_category_description(outfile):
    codes = aaindex1.record_codes()
    indeces = {}
    for code in codes:
        indeces[code] = [aaindex1[code]['category'], aaindex1[code]['description']]

    #save to a csv with columns for code, category, and description
    df = pd.DataFrame(indeces)
    df = df.T
    df.columns = ["category", "description"]
    df.to_csv(outfile)

def produce_supermatrix(lui_list, aln_files, output_file):
    """
    concatenate multiple alignments into one
    :param aln_files: list of fasta files
    :param output_file: output file

    #the goal is to build a supermatrix, which has all the alignments for each sample end-to-end so I can run the gwas on one file 

    #so if we have aln_A 
    >1 
    AAAA
    >2
    AA-A

    and aln_B
    >1
    T--T
    >2
    TTTT

    We will  produce the supermatrix by concatenating the alignments
    super
    >1
    AAAAT--T
    >2
    AA-ATTTT

    #and save a file which maps sites to genes  
    1, aln_A
    2, aln_A
    3, aln_A
    4, aln_A
    5, aln_B
    6, aln_B
    7, aln_B
    8, aln_B
    """

    #initialize the supermatrix
    supermatrix = {}
    for lui in lui_list:
        supermatrix[lui] = ""

    #write the supermatrix out 
    with open("dict.csv", mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=supermatrix.keys())
        writer.writeheader()
        writer.writerows([dict(zip(supermatrix.keys(), values)) for values in zip(*supermatrix.values())])

    #initialize the mapping
    mapping = pd.DataFrame(columns=["site", "gene"])
    #iterate over the alignments
    for aln_file in tqdm(aln_files):
        filtered_records = [
            record
            for record in SeqIO.parse(aln_file, "fasta")
            if record.id.split("|")[0] in lui_list
        ]
        #get the alignment
        for record in filtered_records:
            #get the sample name
            sample_name = record.id.split("|")[0]
            #get the sequence
            seq = str(record.seq)
            #add the sequence to the supermatrix
            supermatrix[sample_name] += seq
        
        #get the gene name
        gene_name = os.path.basename(aln_file).split("_")[0]
        #add the gene name to the mapping
        for i in range(len(seq)):
            #add the site to the mapping
            mapping = pd.concat([mapping, pd.DataFrame({"site": [i], "gene": [gene_name]})], ignore_index=True)
        #pad the supermatrix with gaps if any sequence is shorter than the longest sequence
        max_length = max(len(seq) for seq in supermatrix.values())
        for lui in supermatrix.keys():
            if len(supermatrix[lui]) < max_length:
                #pad with gaps
                supermatrix[lui] += '-' * (max_length - len(supermatrix[lui]))
    #save the supermatrix
    with open(output_file, 'w') as f:
        for lui in lui_list:
            f.write(f">{lui}\n")
            f.write(supermatrix[lui] + "\n")

    #save the mapping
    mapping.to_csv(output_file + "_mapping.csv", index=False)

def get_sites_from_alignment_fasta(aln_file, samples_per_site=0.9):
    """
    prepare an alignment for an association study
    :param aln_file: path to alignment file
    :param samples_per_site: minimum proportion of samples with a site to keep it

    :return: a pandas DataFrame with the alignment, where each column is a site and each row is a sample
    #the column names should be fastaFileBasename_idx 
    """
    # read the alignment
    aln = list(SeqIO.parse(aln_file, 'fasta'))
    # get the alignment length
    aln_len = len(aln[0].seq)

    #assert that all sequences have the same length
    for record in aln:
        assert len(record.seq) == aln_len, 'Sequences have different lengths'

    # create a dictionary to store the alignment
    aln_dict = {i: [''] * aln_len for i in range(len(aln))}

    # fill the dictionary
    for i, record in tqdm(enumerate(aln)):
        aln_dict[i] = list(record.seq.upper())

    # create a DataFrame
    aln_df = pd.DataFrame(aln_dict).T

    # rename the columns
    fastabasename = os.path.basename(aln_file).split(".")[0]
    aln_df.columns = [f"{fastabasename}_{i}" for i in aln_df.columns]

    #rename the rows
    aln_df.index = [record.id.split("_")[0] for record in aln]

    # remove sites with too many missing data (encoded as "-", or gaps )
    print("Removing sites with too many missing data")
    print("Before removing sites:", aln_df.shape)
    #remove sites with too many missing data
    aln_df = aln_df.loc[:, (aln_df == '-').mean() < (1 - samples_per_site)]

    print("After removing sites:", aln_df.shape)

    return aln_df

def encode_using_pcas(aln_df, pcs,pc_limit=3):
    """
    :param aln df: as above, which shall have rows as samples and columns as sites
    :param pca_path: path to the PCA results, which should be a csv file with the first column as the amino acid names and the rest as the PCA components

    we want to expand the alignment df of n samples * m sites into n samples * (m * pca_components), 
    where pca_components is the number of PCA components to use. 

    So if the first three PCS of W are 1.2, 0.8, and 0.5 and for A it is 0.2,0.3,0.1 and we input the alignment:
    W A
    W W

    We will want the output to be:
    1.2 0.8 0.5 0.2 0.3 0.1
    1.2 0.8 0.5 1.2 0.8 0.5

    with sitenames as site_pc1, site_pc2, etc.
    """
    #initialize the encoded df (n samples * (m sites * pca_components))
    encoded_df = pd.DataFrame(index=aln_df.index)
    #get the amino acid names from the pca file
    aa_names = pcs["AminoAcid"].values
    #get the pca components
    pca_components = pcs.iloc[:, 0:pc_limit].values

    def get_pca_components(aa):
        """
        get the pca components for an amino acid
        :param aa: amino acid
        :return: pca components
        """
        allowed_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        if aa not in allowed_aa:
            #if the amino acid is a gap, return a vector of zeros
            return np.zeros(pc_limit)
        #get the index of the amino acid in the pca file
        aa_index = list(aa_names).index(aa)
        #get the pca components for the amino acid
        pca = pca_components[aa_index]
        return pca

    encoded_df = aln_df.map(get_pca_components)
    expanded_df = pd.DataFrame(index=encoded_df.index)

    for col in tqdm(encoded_df.columns):
        #get the pca components for the site
        pcas = encoded_df[col].apply(pd.Series)
        pcas.columns = [f"{col}_pc{i+1}" for i in range(pc_limit)]
        expanded_df = pd.concat([expanded_df, pd.DataFrame(pcas)], axis=1)
        

    """
    #set the index to the same as the input df
    encoded_df.index = aln_df.index
    #iterate over the columns, applying the pca components to each sample and each site
    print("unexpanded,",encoded_df)
    #expand the dataframe,
    #so that each site has its own columns for each pca component
    expanded_df = pd.DataFrame(index=encoded_df.index)
    for col in tqdm(encoded_df.columns):
        print("COL",col)
        #get the pca components for the site
        pcas = encoded_df[col].values
        print("PCAS",pcas)
        #expand the pcas
        #expand the pcas into a dataframe
        print(pcas[0])
        assert pc_limit == len(pcas[0]), f"pc_limit {pc_limit} does not match the number of pca components {len(pcas[0])}"
        for i in range(len(pcas)):
            #convert the pca components to a dataframe
            pcas[i] = pd.Series(pcas[i])


        expanded_df = pd.concat([expanded_df, pd.DataFrame(pcas)], axis='index')
        #rename the columns
        expanded_df.columns = [f"{col}_pc{i+1}" for i in range(pcas.shape[1])]

    #set the index to the same as the input df
    expanded_df.index = encoded_df.index
    print("expanded,",expanded_df)
    """
    return expanded_df

if __name__ == "__main__":

    outfile = "data/aaindex1.csv"
    #download_all_idxs_to_csv(outfile)

    aa_pcs = pd.read_csv("data/aaPCs.csv")

    master_df = pd.read_csv("data/tmp/rbcL_aln/merged_aa_counts.csv")
    lui_list = master_df["ID"].unique().tolist()

    aln_files = glob.glob("data/tmp/alignedGenes/*AA_aligned.fasta")

    produce_supermatrix(lui_list, aln_files, "data/tmp/aa_supermatrix.fasta")

    aln_df = get_sites_from_alignment_fasta("data/tmp/aa_supermatrix.fasta", samples_per_site=0.99)
    #expand the supermatrix

    expanded_df = encode_using_pcas(aln_df, aa_pcs, pc_limit=3)

    #save the expanded df
    expanded_df.to_csv("data/tmp/aa_supermatrix_expanded_3pcs.csv", index=True)

    expanded_df = encode_using_pcas(aln_df, aa_pcs, pc_limit=10)

    #save the expanded df
    expanded_df.to_csv("data/tmp/aa_supermatrix_expanded_10pcs.csv", index=True)
