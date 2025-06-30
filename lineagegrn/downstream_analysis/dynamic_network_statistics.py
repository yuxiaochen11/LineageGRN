import os
import ast
import pandas as pd


#Descriptive statistics of the dyngiamic evolution of gene regulatory networks along fate maps
def get_regulators_for_target_gene(target_gene_id: str, node_id: str, input_path: str, regulator_names: list[str]):
    """
    Retrieve regulator information for a specific target gene.

    Parameters:
    target_gene_id (str): The ID of the target gene.
    node_id (str): The ID of the node.
    input_path (str): The path to the input files.
    regulator_names (list[str]): A list of names of regulators.

    Returns:
    dict: A dictionary containing the regulator information, with keys as regulator names 
          and values as their corresponding lambda values.
    """
    file_path = input_path + '/' +'target_gene_' + target_gene_id + '.csv'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [ast.literal_eval(line.strip()) for line in lines]
        output_df = pd.DataFrame(data)
        output_df.columns = ['target_gene_id', 'node_id', 'lambda'] + regulator_names
        grn = output_df.iloc[:, 3:]
        target_id_node_id_grn = grn.loc[output_df.loc[:, 'node_id'] == node_id, :]
    
    return target_id_node_id_grn.iloc[0].to_dict()

def get_targets_for_regulator_gene(regulator_id: str, node_id: str, input_path: str, regulator_names: list[str]):
    """
    Retrieve target genes for a specific regulator.

    Parameters:
    regulator_id (str): The ID of the regulator.
    node_id (str): The ID of the node.
    input_path (str): The path to the input files.
    regulator_names (list[str]): A list of names of regulators.

    Returns:
    dict: A dictionary where the keys are target gene IDs and the values are 
          the corresponding lambda values for the specified regulator.
    """
    target_gene_dict = {}
    
    for filename in os.listdir(input_path):
        file_path = os.path.join(input_path, filename)
        
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            data = [ast.literal_eval(line.strip()) for line in lines]
            output_df = pd.DataFrame(data)
            output_df.columns = ['target_gene_id', 'node_id', 'lambda'] + regulator_names
            
            grn = output_df.iloc[:, 3:]
            target_gene_id = output_df.loc[0, 'target_gene_id']
            target_id_node_id_grn = grn.loc[output_df.loc[:, 'node_id'] == node_id, :]
            
            target_gene_dict.update({target_gene_id: float(target_id_node_id_grn.loc[:, regulator_id])})
    
    return target_gene_dict
