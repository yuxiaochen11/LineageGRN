import os
import sys
import json
import itertools
import pandas as pd 
import numpy as np
from copy import deepcopy
from typing import Dict, List
from scipy.optimize import minimize
from multiprocessing import Pool, cpu_count


from lineagegrn.utils.constant import * 
from lineagegrn.cell_fate_map import *
from lineagegrn.utils.logs import get_logger
from lineagegrn.utils.basic import ATACSample, ATACData, RNAData, generate_key



CPU_NUM = cpu_count()
logger = get_logger()


def get_grn_initial_data(atac_node_data: 'Dict[str, ATACSample]', reg_genes: 'List[str]', target_gene_id: 'List[str]'):
    """
    Initialize gene regulatory network data using ATAC node data.

    Parameters:
    atac_node_data (Dict[str, ATACSample]): A dictionary mapping keys to ATACSample objects.
    reg_genes (List[str]): List of regulatory gene identifiers.
    target_gene_id (List[str]): List of target gene identifiers.

    Returns:
    np.ndarray: An array of initial GRN values for the given genes.
    """
    data_indexer = {idx: gene_id for idx, gene_id in enumerate(reg_genes)}
    data = np.zeros(shape=len(reg_genes), dtype=np.float32)
    for idx, _ in enumerate(data):
        key = generate_key(data_indexer[idx], target_gene_id)
        data[idx] = atac_node_data.get(key).value if atac_node_data.get(key) else 0.0 

    return data

def compute_probility(a, b, p):
    """
    Compute probability based on the comparison of two values.

    Parameters:
    a (float): First value for comparison.
    b (float): Second value for comparison.
    p (float): Base probability value.

    Returns:
    float: Adjusted probability based on the comparison.
    """
    if abs(b) > 0.1:
        b = 1
    else:
        b = 0
    if a == b:
        return p
    else:
        return 1 - p

def generate_ancestor_state_space(ancestor_nodes: 'List[str]'):
    """
    Generate state space array for ancestor nodes with binary combinations.

    Parameters:
    ancestor_nodes (List[str]): List of ancestor node identifiers.

    Returns:
    np.ndarray: Array representing all possible binary state combinations of ancestor nodes.
    """
    results = np.array([[value for _, value in zip(ancestor_nodes, combination)] 
                        for combination in itertools.product([0, 1], repeat=len(ancestor_nodes))])
    return results

def generate_ancestor_state_space1(ancestor_nodes: 'List[str]'):
    """
    Generate a list of dictionaries representing binary state combinations for ancestor nodes.

    Parameters:
    ancestor_nodes (List[str]): List of ancestor node identifiers.

    Returns:
    List[Dict[str, int]]: List where each element is a dictionary mapping each ancestor node to a binary state.
    """
    results = [{node_id: value for node_id, value in zip(ancestor_nodes, combination)} 
               for combination in itertools.product([0, 1], repeat=len(ancestor_nodes))]
    return results

def _parent_grn_inference(edges: 'List[FateMapEdge]', grn_values: 'Dict[str, List[float]]', regulator_gene_nums: 'int'):
    """
    Infer the parent gene regulatory network values based on child edges and GRN values.

    Parameters:
    edges (List[FateMapEdge]): List of edges from parent to children.
    grn_values (Dict[str, List[float]]): Dictionary containing GRN values for child nodes.
    regulator_gene_nums (int): Number of regulator genes.

    Returns:
    List[float]: Inferred GRN values for the parent node.
    """
    child_weight_1 = edges[0].weight
    child_lambda_1 = grn_values[edges[0].end][0]
    child_grn_value_1 = grn_values[edges[0].end][1:]

    child_weight_2 = edges[1].weight
    child_lambda_2 = grn_values[edges[1].end][0]
    child_grn_value_2 = grn_values[edges[1].end][1:]

    parent_grn_value = []

    for j in range(regulator_gene_nums):
        prob1 = np.exp(-np.array(child_lambda_1) * np.array(child_weight_1))
        prob2 = np.exp(-np.array(child_lambda_2) * np.array(child_weight_2))
        if (child_grn_value_1[j] != child_lambda_2 and child_grn_value_1[j] * child_grn_value_2[j] <= 0):
            if prob1 > prob2:
                parent_grn_value.append(child_grn_value_1[j])
            else:
                parent_grn_value.append(child_grn_value_2[j])
        elif child_grn_value_1[j] * child_grn_value_2[j] > 0:
            parent_grn_value.append((prob1) / (prob1 + prob2) * child_grn_value_1[j] + (prob1) / (prob1 + prob2) * child_grn_value_2[j])   

    parent_lambda = (child_lambda_1 + child_lambda_2) / 2   
    return [parent_lambda] + parent_grn_value

class GRNInference:
    def __init__(self, atac_file_path, expression_file_path, fate_map: 'FateMap', saved_dir: 'str'):
        """
        Initialize GRNInference with file paths, fate map, and directory for saving results.

        Parameters:
        atac_file_path: Path to the ATAC data file.
        expression_file_path: Path to the expression data file.
        fate_map (FateMap): Fate map object containing network structure.
        saved_dir (str): Directory path to save inference results.
        """
        self.fate_map = fate_map
        self.atac_data = ATACData(atac_file_path)
        self.rna_data = RNAData(expression_file_path)
        self.saved_dir = saved_dir

    def _likelihood(self, leaves_grn_data, regulator_data, target_data):
        """
        Calculate the likelihood of GRN parameters given gene expression data.

        Parameters:
        leaves_grn_data: Array of GRN data for leaf nodes.
        regulator_data: Array of regulator gene expression data.
        target_data: Array of target gene expression data.

        Returns:
        float: Calculated likelihood value weighted by a constant.
        """
        G = np.reshape(leaves_grn_data, newshape=[len(self.fate_map.node_leaves), len(self.rna_data.regulator_genes)])
        G = G[..., np.newaxis]
        D = target_data - np.einsum('ijk,ikl -> ijl', regulator_data, G)
        exp_sum = np.sum(-0.5 * SIGMA * np.einsum('ijk,ijk -> i', D, D)) / len(self.rna_data.cell_ids)
        constants_sum = -0.5 * len(self.rna_data.regulator_genes) * np.log(2 * np.pi) * len(self.fate_map.node_leaves) 
        result = exp_sum + constants_sum
        return result * LIKELIHOOD_WEIGHT

    def _atac_prior(self, leaves_grn_data: 'np.ndarray', mdata: 'str'):
        """
        Calculate ATAC prior probability based on GRN and metadata.

        Parameters:
        leaves_grn_data (np.ndarray): Array of GRN data for leaf nodes.
        mdata (str): Processed metadata related to ATAC.

        Returns:
        float: Prior probability weighted by a constant.
        """
        leaves_grn_data = np.stack((abs(leaves_grn_data) > 0.1, abs(leaves_grn_data) <= 0.1), axis=-1, dtype=np.int32)
        G = np.reshape(leaves_grn_data, newshape=[len(self.fate_map.node_leaves), len(self.rna_data.regulator_genes), 2])
        p = np.sum(np.einsum('ijk,ijk -> ij', mdata, G))
        return p * ATAC_PRIOR_WEIGHT

    def _lineage_prior(self, leaves_grn_data: 'np.ndarray', lambda_1: 'float', root_grn_data: 'np.ndarray', ancestor_nodes_space: 'np.ndarray', weights: 'np.ndarray'):
        """
        Calculate lineage prior probability based on GRN data across the lineage tree.

        Parameters:
        leaves_grn_data (np.ndarray): Array of GRN data for leaf nodes.
        lambda_1 (float): Hyperparameter controlling weight adjustments.
        root_grn_data (np.ndarray): GRN data at the root node.
        ancestor_nodes_space (np.ndarray): Possible state space of ancestor nodes.
        weights (np.ndarray): Weights array for lineage prior computation.

        Returns:
        float: Logarithm of the lineage prior probability weighted by a constant.
        """
        lambda_1 = np.clip(lambda_1, LAMBDA_LOW_BOUND, LAMBDA_HIGH_BOUND)
        weights = np.power(weights, lambda_1)
        m = 0
        leaf_node_ids = self.fate_map.node_leaves
        ancestor_node_ids = [self.fate_map.node_root] + self.fate_map.node_internals

        leaves_grn_data = np.reshape(leaves_grn_data, newshape=(len(leaf_node_ids), -1)).T

        root_flags = ancestor_nodes_space[:, ancestor_node_ids.index(self.fate_map.node_root)]
        root_flags = np.concatenate([root_flags[..., np.newaxis], 1 - root_flags[..., np.newaxis]], axis=-1)

        leaf_node_space = np.repeat(leaves_grn_data[:, np.newaxis, ...], len(ancestor_nodes_space), axis=1)
        ancestor_nodes_spaces = np.repeat(ancestor_nodes_space[np.newaxis, ...], len(self.rna_data.regulator_genes), axis=0)
        nodes_space = np.array(np.abs(np.concatenate((ancestor_nodes_spaces, leaf_node_space), axis=-1)) > 0.1)
        nodes_space = np.int32(nodes_space[..., np.newaxis] == nodes_space[..., np.newaxis, :])[..., np.newaxis]
        nodes_space = np.concatenate((nodes_space, 1 - nodes_space), axis=-1)   
        values = np.einsum('rnijk,nijk -> rnij', nodes_space, weights)
        values = np.reshape(values + (values == 0), newshape=(len(self.rna_data.regulator_genes), len(ancestor_nodes_space), -1))
        values = np.prod(values, axis=-1)[..., np.newaxis]
        joint_probility = np.einsum('r,rjk -> rjk', abs(root_grn_data), values)

        root_prior1 = np.concatenate([joint_probility, (1 - joint_probility)], axis=-1)
        p = np.einsum('rjk,jk -> rj', root_prior1, root_flags)
        m = np.sum(np.log(np.sum(p, axis=-1) + 1e-10))
        return m * LINEAGE_PRIOR_WEIGHT

    def _get_strength(self, target_gene_id: 'int' = 0, saved: 'bool' = True):
        """
        Estimate gene regulatory network values for a given target gene.

        Parameters:
        target_gene_id (int): Identifier of the target gene to estimate GRN for.
        saved (bool): Flag indicating whether to save the results.

        Returns:
        Dict[str, List[float]]: Inferred GRN values keyed by node identifiers.
        """
        grn_values = {}
        x = list(self.infer_leaf_grn(target_gene_id))

        leaf_node_ids = self.fate_map.node_leaves
        for idx, node_id in enumerate(leaf_node_ids):
            grn_value = x[1+idx*len(self.rna_data.regulator_genes):1+(idx+1)*len(self.rna_data.regulator_genes)]
            grn_values[node_id] = x[:1] + grn_value

        grn_values = self.infer_ancestor_grn(grn_values, target_gene_id)
        if saved:
            self._save(grn_values, target_gene_id)
            logger.info(f"Saved grn values for target_gene_id:{target_gene_id}")

        return grn_values

    def _reshape_likehood_input(self, target_gene_id):
        """
        Reshape RNA data for likelihood computation.

        Parameters:
        target_gene_id: Identifier for the target gene.

        Returns:
        Tuple[np.ndarray, np.ndarray]: Regulator and target data arrays reshaped for likelihood.
        """
        regulator_data = np.array([self.rna_data.get_values(REGULATOR_GENE, node_id).values for node_id in self.fate_map.node_leaves])
        target_data = np.array([self.rna_data.get_values(TARGET_GENE, node_id=node_id)[[target_gene_id]] for node_id in self.fate_map.node_leaves])
        return regulator_data, target_data

    def _reshape_atac_prior_input(self, target_gene_id):
        """
        Reshape ATAC data for prior computation.

        Parameters:
        target_gene_id: Identifier for the target gene.

        Returns:
        np.ndarray: Processed metadata array for ATAC prior calculation.
        """
        samples = [self.atac_data.get_node_data(node_id=node_id) for node_id in self.fate_map.node_leaves]
        mdata = np.abs([self.atac_data.reshape(sample, self.rna_data.regulator_genes, self.rna_data.target_genes)[[target_gene_id]].values for sample in samples])
        component1 = -np.log(1 + np.exp(-BETA_1 - BETA_2 * mdata))
        component2 = (-BETA_1 - BETA_2 * mdata) - np.log(1 + np.exp(-BETA_1 - BETA_2 * mdata))
        mdata = np.concatenate((component1, component2), axis=-1)
        return mdata

    def _reshape_lineage_prior_input(self):
        """
        Prepare inputs for lineage prior computation by generating ancestor state space and weight matrices.

        Returns:
        Tuple[np.ndarray, np.ndarray]: Ancestor state space array and corresponding weights.
        """
        ancestor_node_ids = [self.fate_map.node_root] + self.fate_map.node_internals
        ancestor_nodes_space = generate_ancestor_state_space(ancestor_node_ids)

        node_ids = [node_id for node_id in self.fate_map.nodes.keys()]

        weights = pd.DataFrame(0.0, index=node_ids, columns=node_ids)
        for grn_edge in self.fate_map.edges:
            weights.loc[grn_edge.start, grn_edge.end] = np.exp(-grn_edge.weight)

        weights = weights.values[..., np.newaxis]
        weights = np.concatenate((weights, 1 - weights), axis=-1)
        weights = np.repeat(weights[np.newaxis, ...], len(ancestor_nodes_space), axis=0)

        return ancestor_nodes_space, weights

    def infer_leaf_grn(self, target_gene_id: 'str' = None):
        """
        Infer GRN values for all leaf nodes for a specific target gene.

        Parameters:
        target_gene_id (str): Identifier of the target gene.

        Returns:
        np.ndarray: Optimized GRN values after inference.
        """
        logger.info(f"Start fitting target_gene_id:{target_gene_id}")
        atac_root_data = self.atac_data.get_node_data(node_id=self.fate_map.node_root)
        reg_genes = self.rna_data.regulator_genes
        root_grn_initial = get_grn_initial_data(atac_root_data, reg_genes, target_gene_id)

        leaf_grn_initial = []
        for node_id in self.fate_map.node_leaves:
            node_atac_data = self.atac_data.get_node_data(node_id)
            intx = get_grn_initial_data(node_atac_data, reg_genes, target_gene_id)
            leaf_grn_initial = leaf_grn_initial + list(intx)

        regulator_data, target_data = self._reshape_likehood_input(target_gene_id)
        mdata = self._reshape_atac_prior_input(target_gene_id)
        ancestor_nodes_space, weights = self._reshape_lineage_prior_input()

        def _loss(x):
            y1 = self._lineage_prior(x[1:], x[0], root_grn_initial, ancestor_nodes_space, weights)
            y2 = self._atac_prior(x[1:], mdata)
            y3 = self._likelihood(x[1:], regulator_data, target_data)
            y = (-y1 - y2 - y3)
            return y

        x0 = np.array([LAMBDA_1] + leaf_grn_initial)
        bounds = [(LAMBDA_LOW_BOUND, LAMBDA_HIGH_BOUND)] + [(REGULATION_STRENGTH_LOW_BOUND, REGULATION_STRENGTH_HIGH_BOUND)] * len(leaf_grn_initial)
                                                
        result = minimize(fun=_loss,
                          x0=x0,
                          method=OPTIMIZE_METHOD,
                          bounds=bounds,
                          )

        logger.info(f"Finish inferencing leaves grn value for target_gene_id:{target_gene_id}")
        return result.x

    def infer_ancestor_grn(self, grn_values: 'Dict[str, List[float]]', target_gene_id: 'str'):
        """
        Infer ancestor GRN values by aggregating information from leaf nodes upward.

        Parameters:
        grn_values (Dict[str, List[float]]): Dictionary of GRN values for leaf nodes.
        target_gene_id (str): Identifier of the target gene.

        Returns:
        Dict[str, List[float]]: Updated dictionary including ancestor GRN values.
        """
        stacks = deepcopy(self.fate_map.node_leaves)

        def _isin_stacks(node_ids, stacks):
            flag = False if sum([0 if node_id in stacks else 1 for node_id in node_ids]) else True
            return flag

        idx = 0
        while idx < len(stacks):
            node_id = stacks[idx]
            parent_node_id = self.fate_map.nodes[node_id].upstream_node_id
            if not parent_node_id or parent_node_id in stacks:
                idx += 1
                continue

            parent_node = self.fate_map.get_parent_node(node_id)
            edges = parent_node.directed_edges
            children_ids = [edge.end for edge in edges]

            if not _isin_stacks(children_ids, stacks):
                idx += 1
                continue

            parent_grn_value = _parent_grn_inference(edges, grn_values, len(self.rna_data.regulator_genes))
            stacks.append(parent_node.node_id)
            grn_values[parent_node.node_id] = parent_grn_value
            idx += 1

        logger.info(f"Finish inferencing leaves grn value for target_gene_id:{target_gene_id}")
        return grn_values

    def _save(self, grn_values: 'Dict[str, List[float]]', target_gene_id: 'str'):
        """
        Save the inferred GRN values to a CSV file.

        Parameters:
        grn_values (Dict[str, List[float]]): Dictionary of GRN values keyed by node.
        target_gene_id (str): Identifier of the target gene.
        """
        if not os.path.exists(self.saved_dir):
            os.makedirs(self.saved_dir)

        saved_path = os.path.join(self.saved_dir, f"target_gene_{target_gene_id}.csv")
        with open(saved_path, 'w+') as f:
            for node_id, grn_value in grn_values.items():
                line = [target_gene_id, node_id] + grn_value
                line = json.dumps(line)
                f.write(line + '\n')

    def infer_grn(self, max_processing=CPU_NUM):
        """
        Estimate GRN values for all target genes using parallel processing if specified.

        Parameters:
        max_processing: Maximum number of parallel processes to use.
        """
        if max_processing > 1:
            with Pool(max_processing) as p:
                p.map(self._get_strength, self.rna_data.target_genes)
        else:
            for target_gene_id in self.rna_data.target_genes:
                self._get_strength(target_gene_id)

    def get_target_networks(self, threshold: 'float', reverse=True) -> 'Dict[str, Dict[str, Dict]]':
        """
        Retrieve gene regulatory networks for each target gene based on a threshold.

        Parameters:
        threshold (float): Threshold for filtering GRN values.
        reverse (bool): Flag to reverse direction if needed (unused in logic).

        Returns:
        Dict[str, Dict[str, Dict]]: Nested dictionary of GRN networks keyed by node and target gene.
        """
        grn_dict = {}
        filenames = [f'target_gene_{target_gene_id}.csv' for target_gene_id in self.rna_data.target_genes]
        for filename in filenames:
            file_path = os.path.join(self.saved_dir, filename)
            with open(file_path, 'r') as f:
                for line in f.readlines():
                    line = json.loads(line.strip())
                    target_gene_id = line[0]
                    node_id = line[1]
                    lambda_1 = line[2]
                    grn_value = [0.0 if abs(value) < threshold else value for value in line[3:]]
                    if node_id not in grn_dict:
                        grn_dict[node_id] = {target_gene_id: {'lambda': lambda_1, 'grn_value': grn_value}}
                    else:
                        grn_dict[node_id][target_gene_id] = {'lambda': lambda_1, 'grn_value': grn_value}

        return grn_dict
