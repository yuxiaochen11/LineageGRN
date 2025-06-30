import os
import re
import csv
import pandas as pd
import numpy as np

def run_muscle(input_fasta, output_aln):
    """
    Perform multiple sequence alignment using MUSCLE.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_aln (str): Path to the output alignment file (default: muscle_output.aln).

    Returns:
        str: Path to the MUSCLE output alignment file.
    """
    try:
        command = f"muscle -in {input_fasta} -out {output_aln}"
        os.system(command)
        print(f"MUSCLE alignment complete. Results saved to: {output_aln}")
        return output_aln
    except Exception as e:
        print("MUSCLE execution failed:", e)
        return None

def parse_muscle_output(aligned_file):
    """
    Parse the MUSCLE alignment output file.

    Parameters:
        aligned_file (str): Path to the MUSCLE output alignment file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are aligned sequences.
    """
    if aligned_file is None:
        print("No output file to parse.")
        return None

    try:
        aligned_sequences = {}
        with open(aligned_file, 'r') as f:
            seq_id = None
            for line in f:
                if line.startswith('>'):
                    seq_id = line.strip()[1:]
                    aligned_sequences[seq_id] = ''
                elif seq_id:
                    aligned_sequences[seq_id] += line.strip()

        return aligned_sequences
    except Exception as e:
        print(f"Error parsing alignment file: {e}")
        return None

def save_aligned_sequences(aligned_sequences, output_file):
    """
    Save parsed aligned sequences to a CSV file.

    Parameters:
        aligned_sequences (dict): Dictionary of sequence IDs and aligned sequences.
        output_file (str): Path to the output CSV file.
    """
    try:
        df = pd.DataFrame.from_dict(aligned_sequences, orient='index', columns=['Aligned Sequence'])
        df.index.name = 'Sequence ID'
        df.to_csv(output_file, index=True)
        print(f"Aligned sequences saved to: {output_file}")
    except Exception as e:
        print(f"Error saving aligned sequences: {e}")

def csv_to_fasta(input_csv, output_fasta):
    """
    Convert a CSV file containing sequence data to FASTA format.

    Parameters:
        input_csv (str): Path to the input CSV file.
        output_fasta (str): Path to the output FASTA file.
    """
    try:
        with open(input_csv, 'r') as csvfile, open(output_fasta, 'w') as fastafile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                cell_id = row['cell_id']
                raw_barcode = row['raw_barcodes']
                fastafile.write(f">{cell_id}\n{raw_barcode}\n")
        print(f"FASTA file generated: {output_fasta}")
    except Exception as e:
        print(f"Error converting CSV to FASTA: {e}")

def sort_by_sequence_id(input_file, output_file):
    """
    Sort a CSV file by the 'Sequence ID' column and save the sorted file.

    Parameters:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output sorted CSV file.
    """
    try:
        df = pd.read_csv(input_file)
        if 'Sequence ID' not in df.columns:
            raise ValueError("The file lacks a 'Sequence ID' column.")

        df['ID_Integer'] = df['Sequence ID'].apply(lambda x: int(re.search(r'_(\d+)$', x).group(1)))
        df_sorted = df.sort_values(by='ID_Integer').drop(columns=['ID_Integer'])
        df_sorted.to_csv(output_file, index=False)
        print(f"File sorted by 'Sequence ID' and saved to: {output_file}")
    except Exception as e:
        print(f"Error sorting file: {e}")

def read_aligned_sequences(file_path):
    """
    Read aligned sequences from a CSV file.

    Parameters:
        file_path (str): Path to the CSV file containing aligned sequences.

    Returns:
        DataFrame: A pandas DataFrame containing sequence IDs and aligned sequences.
    """
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def generate_mutation_matrix(df):
    """
    Generate a mutation matrix based on sequence differences in specified regions.

    Parameters:
        df (DataFrame): DataFrame containing sequence IDs and aligned sequences.

    Returns:
        DataFrame: Mutation matrix with regions as columns and mutation types as values.
    """
    sequences = df['Aligned Sequence']
    sequence_ids = df['Sequence ID']

    window_lengths = [38, 18, 14, 12, 15, 7, 9, 7, 10, 14, 22]
    seq_length = len(sequences.iloc[0])

    if sum(window_lengths) != seq_length:
        raise ValueError("Window lengths do not match sequence length.")

    reference_sequence = sequences.iloc[0]
    ref_cut = [reference_sequence[start:start + length] for start, length in zip(np.cumsum([0] + window_lengths[:-1]), window_lengths)]

    region_mutation_dicts = [{} for _ in range(len(window_lengths))]
    mutation_matrix = []

    for seq in sequences:
        mutations = []
        for idx, region_seq in enumerate([seq[start:start + length] for start, length in zip(np.cumsum([0] + window_lengths[:-1]), window_lengths)]):
            if region_seq == ref_cut[idx]:
                mutations.append(0)
            else:
                mutation_dict = region_mutation_dicts[idx]
                if region_seq in mutation_dict:
                    mutations.append(mutation_dict[region_seq])
                else:
                    mutation_id = len(mutation_dict) + 1
                    mutation_dict[region_seq] = mutation_id
                    mutations.append(mutation_id)
        mutation_matrix.append(mutations)

    region_labels = [f"Region_{i+1}" for i in range(len(window_lengths))]
    mutation_df = pd.DataFrame(mutation_matrix, columns=region_labels)
    mutation_df.insert(0, 'Sequence ID', sequence_ids)
    return mutation_df

def save_mutation_matrix(mutation_df, output_file):
    """
    Save the mutation matrix to a CSV file.

    Parameters:
        mutation_df (DataFrame): DataFrame containing the mutation matrix.
        output_file (str): Path to the output CSV file.
    """
    try:
        mutation_df.to_csv(output_file, index=False)
        print(f"Mutation matrix saved to: {output_file}")
    except Exception as e:
        print(f"Error saving mutation matrix: {e}")



