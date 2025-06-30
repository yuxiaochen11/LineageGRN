import os
import ast
import numpy as np
import pandas as pd
from io import StringIO
from Bio import Phylo
from tqdm import tqdm
from pandas import DataFrame
from typing import Dict, List, Tuple
from .constant import * 
from .logs import get_logger
logger = get_logger()

def generate_key(regulator_gene, target_gene):
    return f"{regulator_gene}->{target_gene}"

class ATACSample:
    def __init__(self, regulator_gene, target_gene, value, node_id:str = None):
        self.regulator_gene = regulator_gene
        self.target_gene = target_gene
        self.value = value
        self.node_id = node_id
        self.sample_id = generate_key(regulator_gene, target_gene)
    
    def __repr__(self) -> str:
        return f"{self.node_id}\t{self.regulator_gene}\t{self.target_gene}\t{self.value}"

class ATACData:
    def __init__(self, atac_file_path):
        self.data = self._serialize(atac_file_path)
    
    def reshape(self,samples: Dict[str,ATACSample], regulator_gene_ids, target_gene_ids):
        indexes = regulator_gene_ids
        columns = target_gene_ids
        data = pd.DataFrame(0.0, index=indexes, columns=columns)
        for sample in samples.values():
            data.loc[sample.regulator_gene, sample.target_gene] = sample.value

        return data

    def _serialize(self, atac_file_path):
        data = {}
        with open(atac_file_path,'r') as f:
            for line in f.readlines():
                regulator_gene_id, target_gene_id, value, node_id = line.strip().split(',')
                if node_id not in data:
                    data[node_id] = {}
                sample = ATACSample(regulator_gene=regulator_gene_id, target_gene=target_gene_id, value=np.float32(value), node_id=node_id)
                data[node_id][sample.sample_id] = sample
        return data

    def get_sample(self, node_id:str, regulator_gene, target_gene):
        key = generate_key(regulator_gene, target_gene)
        return self.data[node_id][key]
    
    def get_node_data(self, node_id: str)->Dict[str,ATACSample]:
        return self.data[node_id]
    
    @property
    def all_node_ids(self,):
        return list(self.data.keys())
    
    def __repr__(self) -> str:
        return f"Nodes:{self.all_node_ids}"

class GeneSample:
    def __init__(self, gene_id:str, cell_id:str, value:float, gene_type:str, node_id:str):
        self.gene_id = gene_id
        self.cell_id = cell_id
        self.value = value
        self.gene_type = gene_type
        self.node_id = node_id
    
    def __repr__(self) -> str:
        return f"{self.gene_id}\t{self.cell_id}\t{self.value}"
    
    def to_list(self):
        return [self.gene_id, self.cell_id, self.value]


class NodeSample:
    def __init__(self,regulator_gene_data:pd.DataFrame,target_gene_data:pd.DataFrame,node_id:str):
        self._regulator_gene_data = regulator_gene_data
        self._target_gene_data = target_gene_data
        self.node_id = node_id

    @property
    def regulator_data(self):
        return self._regulator_gene_data
    @property
    def target_data(self):
        return self._target_gene_data
    @property
    def cell_ids(self):
        return list(self._regulator_gene_data.index)
    @property
    def regulator_gene_ids(self):
        return list(self._regulator_gene_data.columns)
    @property
    def target_gene_ids(self):
        return list(self._target_gene_data.columns)

class RNAData:
    def __init__(self,expression_file_path) -> None:
        self.data, self.node_ids = self._serialize(expression_file_path)

    def _reshape(self, gene_samples:List[GeneSample], indexes:List[str] = None, columns:List[str] = None)->DataFrame:

        if not indexes:
            indexes = list(set([sample.cell_id for sample in gene_samples]))
        if not columns:
            columns = list(set([sample.gene_id for sample in gene_samples]))

        if not isinstance(indexes, List):
            indexes = list(indexes)

        if not isinstance(columns, List):
            columns = list(columns)

        data = pd.DataFrame(0.0, index=indexes, columns=columns)
        for sample in gene_samples:
            data.loc[sample.cell_id,sample.gene_id] = sample.value

        for column in data.columns:
            values = np.array(data[column])
            _mean = np.mean(values)
            _std = np.std(values)
            data[column] = (values - _mean)/_std if _std else (values - _mean)

        return data.fillna(0.0)

    def _serialize(self,expression_file_path)->Tuple[Dict[str,NodeSample],List[str]]:
        data = {}
        node_ids = []

        gene_ids = {REGULATOR_GENE:[],TARGET_GENE:[]}
        cell_ids = []

        with open(expression_file_path,'r') as f:
            for line in f.readlines():
                gene_id, cell_id, value, node_id, gene_type = line.strip().split(',')
                if gene_id not in gene_ids[gene_type]:
                    gene_ids[gene_type].append(gene_id)
                if cell_id not in cell_ids:
                    cell_ids.append(cell_id)

                gene_sample = GeneSample(gene_id,cell_id,float(value),gene_type,node_id)
                if node_id not in data:
                    data[node_id] = [gene_sample]
                else:
                    data[node_id] += [gene_sample]


        for node_id, gene_samples in data.items(): 
            regulator_gene_samples = [gene_sample for gene_sample in gene_samples if gene_sample.gene_type == REGULATOR_GENE ]
            target_gene_samples = [gene_sample for gene_sample in gene_samples if gene_sample.gene_type == TARGET_GENE ]
            logger.info(f"Serialize node_id {node_id} expression data of {REGULATOR_GENE}")
            regulator_gene_data = self._reshape(regulator_gene_samples,indexes=cell_ids, columns=gene_ids[REGULATOR_GENE])
            logger.info(f"Serialize node_id {node_id} expression data of {TARGET_GENE}")
            target_gene_data = self._reshape(target_gene_samples, indexes=cell_ids, columns=gene_ids[TARGET_GENE])
            node_sample = NodeSample(
                regulator_gene_data=regulator_gene_data,
                target_gene_data=target_gene_data,
                node_id=node_id,
            )
            data[node_id] = node_sample
            node_ids.append(node_id)
        self._gene_ids = gene_ids
        self._cell_ids = cell_ids
        return data, node_ids

    
    @property
    def regulator_genes(self)->List[str]:
        return list(self._gene_ids[REGULATOR_GENE])

    @property
    def target_genes(self)->List[str]:
        return list(self._gene_ids[TARGET_GENE])
    
    @property
    def cell_ids(self):
        return list(self._cell_ids)
    
    def get_values(self, sample_type, node_id:str):
        if sample_type == REGULATOR_GENE:
            data = self.data[node_id].regulator_data        
        elif sample_type == TARGET_GENE:
            data = self.data[node_id].target_data
        else:
            raise ValueError(f"No Sample type {sample_type}")

        return data


def newick_to_edge_length_dict(newick_str):
    """
    Convert a Newick string into a dictionary mapping parent-child node pairs to branch lengths.
    Requires that all internal nodes in the Newick string are explicitly named (e.g., P0, P1...).

    Parameters:
    newick_str (str): A Newick string that includes both branch lengths and internal node names.

    Returns:
    dict: A dictionary with keys in the format 'Parent->Child' and values as branch lengths (float).
    """
    handle = StringIO(newick_str)
    tree = Phylo.read(handle, "newick")
    edge_dict = {}

    def traverse(clade):
        for subclade in clade.clades:
            parent = clade.name
            child = subclade.name
            edge_dict[f"{parent}->{child}"] = subclade.branch_length
            traverse(subclade)

    traverse(tree.root)
    return edge_dict



def get_dynamic_networks(input_path, fate_map, threshold, regulator_names, target_gene_names):
    """
    Retrieve the Gene Regulatory Network (GRN) for each node in the fate map.

    Parameters:
    input_path (str): The path to the input files containing regulatory data.
    fate_map (object): An object that contains node information.
    threshold (float): A threshold value below which regulatory values are set to zero.
    regulator_names (list[str]): A list of names of regulators.
    target_gene_names (list[str]): A list of names of target genes.

    Returns:
    dict: A dictionary where keys are node IDs and values are DataFrames containing 
        the regulatory values for the target genes associated with each node.
    """
    grns_dict = {}
    nodes = fate_map.nodes  # Retrieve the list of nodes from the fate map

    for node in nodes:
        grn_df = pd.DataFrame()  # Initialize an empty DataFrame for the current node
        
        for filename in ['target_gene_'+str(i)+'.csv' for i in target_gene_names]:
            file_path = os.path.join(input_path, filename)
            
            if os.path.isfile(file_path):
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                
                # Read data and create a DataFrame for target gene regulatory data
                data = [ast.literal_eval(line.strip()) for line in lines]
                target_gene_data = pd.DataFrame(data).iloc[:, [0]+list(range(3, pd.DataFrame(data).shape[1]))]
                target_gene_data.columns = ['target_gene_id']+regulator_names
                target_gene_data.index = list(pd.DataFrame(data).iloc[:, 1])
                # Append regulatory data for the current node to the GRN DataFrame
                grn_df = pd.concat([grn_df, pd.DataFrame(target_gene_data.loc[node, :]).T], axis=0)
        grn_df.set_index(grn_df.columns[0], drop=True, inplace=True)
        grn_df = grn_df.applymap(lambda x: 0 if abs(x) < threshold else x)  # Apply threshold
        grn_df = grn_df.loc[target_gene_names,:]
        grns_dict.update({node: grn_df})  # Update the dictionary with the current node's GRN
    return grns_dict


def get_gene_interaction(grns):
    stacked_df = grns.stack().reset_index()
    stacked_df.columns = ['Target', 'Regulator', 'Value']
    stacked_df = stacked_df[['Regulator', 'Target', 'Value']]
    gene_interaction = stacked_df[stacked_df['Value'] > 0].reset_index(drop=True)
    return gene_interaction

def convert_grn_format(grn_index, grn_name, regulatory_genes_name, target_genes_name, output_path):
    grn_index = pd.DataFrame(grn_index)
    grn_index.columns = regulatory_genes_name
    grn_index.index = target_genes_name
    grn_melted = grn_index.reset_index().melt(id_vars='index', var_name='regulatory_gene', value_name='regulation_strength')
    grn_melted.columns = ['target_gene', 'regulatory_gene', 'regulation_strength']
    grn_melted = grn_melted[abs(grn_melted['regulation_strength']) >= 0.05]
    

    #grn_melted['regulation_strength'] = grn_melted['regulation_strength'].apply(lambda a: np.exp(5*a))
    
    grn_melted = grn_melted[grn_melted['target_gene'] != grn_melted['regulatory_gene']]
    grn_melted = grn_melted[['regulatory_gene', 'target_gene', 'regulation_strength']]
    grn_melted.to_csv(output_path + str(grn_name) + '_regulatory_relationships.txt', sep='\t', index=False)
    
    return grn_melted
