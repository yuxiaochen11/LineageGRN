import re
import dendropy
import numpy as np
import pandas as pd
from dendropy import Tree
from typing import Dict, List
from itertools import combinations
from scipy.optimize import minimize_scalar
from scipy.cluster.hierarchy import ClusterNode,linkage,to_tree 
from lineagegrn.utils.logs import get_logger
from lineagegrn.utils.constant import NODE_INTERNAL, NODE_LEAF, NODE_ROOT
logger = get_logger()




class FateMapEdge:
    def __init__(self,start:str, end:str, weight:float) -> None:
        self.start = start
        self.end = end 
        self.weight = weight
    def __repr__(self):
        return f"{self.start}\t{self.end}\t{self.weight}"

class FateMapNode:
    def __init__(self,node_id:str, directed_edges:List[FateMapEdge] = None, upstream_node_id:str = None, node_type:str = None,high_expression_genes_in_leaf:List[str]=None) -> None:
        self.node_id = node_id
        self.node_type = node_type
        self.directed_edges = directed_edges
        self.upstream_node_id = upstream_node_id
        self.high_expression_genes_in_leaf=high_expression_genes_in_leaf

    def __repr__(self):
        return f"{self.node_id}\t{self.node_type}\tDirected_Edges:{len(self.directed_edges)}"

class FateMap:
    def __init__(self,edges:List[FateMapEdge],high_expression_genes:Dict[str,List[str]] = None) -> None:
        self.edges = edges
        self.nodes = self._get_nodes(high_expression_genes)
    
    def _get_nodes(self,high_expression_genes:Dict[str,List[str]])->Dict[str,FateMapNode]:
        nodes:Dict[str,FateMapNode] = {}
        for edge in self.edges:
            start = edge.start
            end = edge.end
            if start not in nodes:
                directed_edges = [edge]
                nodes[start] = FateMapNode(
                    node_id=start,
                    directed_edges=directed_edges,
                    node_type=NODE_ROOT,
                    high_expression_genes_in_leaf=high_expression_genes[end] if high_expression_genes is not None else None
                )
            else:
                directed_edges = nodes[start].directed_edges
                directed_edges.append(edge)
            
            if nodes[start].upstream_node_id:
               nodes[start].node_type = NODE_INTERNAL 

            if end not in nodes:
                directed_edges = []
                nodes[end] = FateMapNode(
                    node_id=end,
                    upstream_node_id=nodes[start].node_id,
                    directed_edges=directed_edges,
                    node_type=NODE_LEAF,
                    high_expression_genes_in_leaf=high_expression_genes[end] if high_expression_genes is not None else None 
                )
            else:
                nodes[end].upstream_node_id = nodes[start].node_id
                

            if nodes[end].directed_edges:
               nodes[end].node_type = NODE_INTERNAL
               
        return nodes

    def get_parent_node(self,child_node_id:str)->FateMapNode:
        parent_node_id = self.nodes[child_node_id].upstream_node_id
        return self.nodes[parent_node_id]    
    
    def get_path(self, leaf_node_id: str) -> List[FateMapNode]:
        path = []
        child = leaf_node_id
        while child is not None:
            path.append(self.nodes[child].node_id)
            child = self.nodes[child].upstream_node_id
        path.reverse()  
        return path

    @property
    def node_leaves(self):
        return [node_id for node_id in self.nodes if self.nodes[node_id].node_type == NODE_LEAF]

    @property
    def node_internals(self):
        return [node_id for node_id in self.nodes if self.nodes[node_id].node_type == NODE_INTERNAL]

    @property
    def node_root(self):
        return [node_id for node_id in self.nodes if self.nodes[node_id].node_type == NODE_ROOT][0]
    
    def __repr__(self):
        return f"Nodes:{len(self.nodes)}, Edges:{len(self.edges)}"

def generate_newick(fate_map:FateMap):
        
        root=fate_map.node_root
        leaves={}
        internal_nodes={}
        edges={root:0}
        
        for node_id in [fate_map.node_root]+fate_map.node_internals+fate_map.node_leaves:
            directed_edge=fate_map.nodes[node_id].directed_edges
            if len(directed_edge)!=0:
                internal_nodes.update({node_id:[directed_edge[0].end,directed_edge[1].end]})
                if directed_edge[0].end in fate_map.node_leaves:
                    leaves.update({directed_edge[0].end:directed_edge[0].weight})
                else:
                    edges.update({directed_edge[0].end:directed_edge[0].weight})
                if directed_edge[1].end in fate_map.node_leaves:
                    leaves.update({directed_edge[1].end:directed_edge[1].weight})
                else:
                    edges.update({directed_edge[1].end:directed_edge[1].weight})
            else:
                next
        
        edges.update(leaves)
    
        def build_subtree(node):
            if node in leaves:  
                return f"{node}:{edges[node]}"
            elif node in internal_nodes: 
                subtrees = [build_subtree(child) for child in internal_nodes[node]]
                return f"({','.join(subtrees)}){node}:{edges.get(node, 0)}"  
            else:
                raise ValueError(f"unknown node: {node}")
        return f"{build_subtree(root)};"



def compute_default_allele_prob_new(sc_mat):
    """
    Generate uniform allele emergence probabilities for each site (column) in a single-cell matrix.

    Parameters
    ----------
    sc_mat : np.ndarray
        A 2D array (cells x sites). Each entry indicates the mutation state for a specific site in a specific cell.

    Returns
    -------
    np.ndarray
        A 2D array of the same shape as `sc_mat`. Each element represents a uniform probability 
        of emergence for the corresponding site.
    """
    obs_allele_probs_matrix = np.zeros((sc_mat.shape[0], sc_mat.shape[1]))

    # For each column (mutation site):
    for j in range(obs_allele_probs_matrix.shape[1]):
        # Extract unique nonzero mutation states in this column
        muts = np.sort(np.unique(sc_mat[:, j]))
        muts = muts[muts != 0]
        if len(muts) == 0:
            # If no nonzero mutations exist, set the probability to 0
            obs_prob = 0.0
        else:
            # Uniform probability among all observed (nonzero) mutations
            obs_prob = 1.0 / len(muts)

        # Fill the column with the uniform probability
        obs_allele_probs_matrix[:, j] = obs_prob * np.ones(sc_mat.shape[0])

    return obs_allele_probs_matrix


def estimate_allele_prob_new(sc_mat):
    """
    Estimate the allele emergence probabilities for each site based on observed frequencies in the single-cell matrix.

    Parameters
    ----------
    sc_mat : np.ndarray
        A 2D array (cells x sites). Each entry indicates the mutation state for a specific site in a specific cell.

    Returns
    -------
    np.ndarray
        A 2D array of the same shape as `sc_mat`. Each element represents the estimated probability 
        that a cell carries a particular nonzero mutation (conditional on the site being mutated).
    """
    obs_allele_probs_matrix = np.zeros((sc_mat.shape[0], sc_mat.shape[1]))

    # For each column (mutation site):
    for j in range(obs_allele_probs_matrix.shape[1]):
        # Extract unique nonzero mutation states
        muts = np.sort(np.unique(sc_mat[:, j]))
        muts = muts[muts != 0]

        # Compute the frequency-based probability for each mutation state
        for mut_val in muts:
            # Probability that a cell has the mutation 'mut_val' among all nonzero states
            prob = np.sum(sc_mat[:, j] == mut_val) / np.sum(sc_mat[:, j] != 0)

            # Update those positions to this probability
            index = (sc_mat[:, j] == mut_val)
            obs_allele_probs_matrix[index, j] = prob

    return obs_allele_probs_matrix


def compute_MRCA_time(t_S, sc_mat, mut_prob, alpha_vector, beta):
    """
    Compute the Most Recent Common Ancestor (MRCA) times for all cell pairs.

    Parameters
    ----------
    t_S : float
        The sampling time.
    sc_mat : np.ndarray
        A 2D array (cells x sites). Each entry indicates the mutation state for a specific site in a specific cell.
    mut_prob : np.ndarray
        A 2D array (cells x sites). Each entry indicates the probability of a mutation at a site for a particular cell.
    alpha_vector : np.ndarray
        A 1D array of alpha parameters (one per site).
    beta : float
        A parameter for the (1 + beta * t)^(-alpha) model (scale factor in the integrals).

    Returns
    -------
    np.ndarray
        A 2D array (cells x cells), where entry (i, j) is the estimated MRCA time between cell i and j.
    """
    def compute_integrate(val, t_S, gene1, gene2, case1, case2, case3, case4, case5, alpha_vector, beta):
        """
        Inner function to compute the negative log-likelihood given a candidate t_MRCA value.
        """
        t_MRCA = val

        # Compute terms (1 + beta*t)^{-alpha_i} for relevant timepoints
        # Note: I1 is always 1, so we don't store it
        I2 = (1 + beta * t_MRCA) ** (-alpha_vector)         # (1 + beta*t_MRCA)^(-alpha)
        I3 = (1 + beta * t_S) ** (-alpha_vector)            # (1 + beta*t_S)^(-alpha)
        I4 = (1 + 2*beta*t_S - beta*t_MRCA) ** (-alpha_vector)

        # case1: Both cells share the same nonzero mutation
        prob_case1 = gene1[case1] + (-gene1[case1] + gene1[case1]**2) * I2[case1] \
                     - 2*(gene1[case1]**2)*I3[case1] + (gene1[case1]**2)*I4[case1]

        # case2: Both cells have no mutation (0)
        prob_case2 = I4[case2]

        # case3: Cells have different nonzero mutations
        prob_case3 = (I2[case3] - 2*I3[case3] + I4[case3]) * gene1[case3] * gene2[case3]

        # case4: First cell mutated, second cell not mutated
        prob_case4 = (I3[case4] - I4[case4]) * gene1[case4]

        # case5: First cell not mutated, second cell mutated
        prob_case5 = (I3[case5] - I4[case5]) * gene2[case5]

        # Sum the log probabilities across all sites
        loglikelihood = (np.sum(np.log(prob_case1)) +
                         np.sum(np.log(prob_case2)) +
                         np.sum(np.log(prob_case3)) +
                         np.sum(np.log(prob_case4)) +
                         np.sum(np.log(prob_case5)))

        # We return negative log-likelihood for minimization
        return -loglikelihood

    MRCA_matrix = np.zeros((sc_mat.shape[0], sc_mat.shape[0]))

    # For each pair of cells (k, j):
    for k in range(sc_mat.shape[0]):
        for j in range(k, sc_mat.shape[0]):
            cell1 = sc_mat[k, :]
            cell2 = sc_mat[j, :]

            gene1 = mut_prob[k, :]
            gene2 = mut_prob[j, :]

            # Define cases for each site
            case1 = ((cell1 == cell2) & (cell1 != 0))
            case2 = ((cell1 == cell2) & (cell1 == 0))
            case3 = ((cell1 != cell2) & (cell1 != 0) & (cell2 != 0))
            case4 = ((cell1 != 0) & (cell2 == 0))
            case5 = ((cell1 == 0) & (cell2 != 0))

            # Minimize negative log-likelihood wrt t_MRCA
            pred_t_M = minimize_scalar(
                compute_integrate,
                args=(t_S, gene1, gene2, case1, case2, case3, case4, case5, alpha_vector, beta),
                method='bounded', bounds=(0.5, 15)
            )

            # Update symmetric matrix
            MRCA_matrix[k, j] = pred_t_M.x
            MRCA_matrix[j, k] = pred_t_M.x

    return MRCA_matrix


def compute_alpha(t_S, sc_mat, mut_prob, beta, MRCA_matrix):
    """
    Compute the alpha parameters for each mutation site, given the MRCA matrix and other model parameters.

    Parameters
    ----------
    t_S : float
        The sampling time.
    sc_mat : np.ndarray
        A 2D array (cells x sites). Each entry indicates the mutation state for a specific site in a specific cell.
    mut_prob : np.ndarray
        A 2D array (cells x sites). Each entry indicates the probability of a mutation at a site for a particular cell.
    beta : float
        A parameter for the (1 + beta * t)^(-alpha) model (scale factor in the integrals).
    MRCA_matrix : np.ndarray
        A 2D array (cells x cells). The MRCA times between each pair of cells.

    Returns
    -------
    np.ndarray
        A 1D array of optimized alpha parameters (one per mutation site).
    """
    def compute_alpha_likelihood(val, t_S, sc_mat, mut_prob, MRCA_matrix, beta, i):
        """
        Inner function to compute the negative log-likelihood for a single alpha_i value, 
        given site i and the MRCA matrix for all cell pairs.
        """
        alpha_i = val
        loglikelihood = 0

        # For each pair of cells (k, j):
        for k in range(sc_mat.shape[0]):
            for j in range(k, sc_mat.shape[0]):
                cell1 = sc_mat[k, i]
                cell2 = sc_mat[j, i]
                gene1 = mut_prob[k, i]
                gene2 = mut_prob[j, i]

                I2 = (1 + beta * MRCA_matrix[k, j]) ** (-alpha_i)
                I3 = (1 + beta * t_S) ** (-alpha_i)
                I4 = (1 + 2*beta*t_S - beta*MRCA_matrix[k, j]) ** (-alpha_i)

                # Evaluate the probability based on the same cases used in compute_MRCA_time
                if cell1 == cell2:
                    if cell1 != 0:
                        prob = (gene1
                                + (-gene1 + gene1**2) * I2
                                - 2 * (gene1**2) * I3
                                + (gene1**2) * I4)
                    else:
                        prob = I4
                else:
                    if (cell1 != 0) and (cell2 != 0):
                        prob = (I2 - 2*I3 + I4) * gene1 * gene2
                    elif (cell1 != 0) and (cell2 == 0):
                        prob = (I3 - I4) * gene1
                    elif (cell1 == 0) and (cell2 != 0):
                        prob = (I3 - I4) * gene2

                loglikelihood += np.log(prob)

        # Return negative log-likelihood
        return -loglikelihood

    alpha_vector = np.zeros(sc_mat.shape[1])

    # Optimize alpha_i for each site i
    for i in range(len(alpha_vector)):
        pred_alpha = minimize_scalar(
            compute_alpha_likelihood,
            args=(t_S, sc_mat, mut_prob, MRCA_matrix, beta, i),
            method='bounded', bounds=(0, 10)  # 10 is an upper bound for the search
        )
        alpha_vector[i] = pred_alpha.x

    return alpha_vector


def alternating_optimization(t_S, sc_mat, mut_prob, beta, max_iter=100, tol=1e-6):
    """
    Perform an alternating optimization of MRCA_matrix and alpha_vector.

    Parameters
    ----------
    t_S : float
        Sampling time.
    sc_mat : np.ndarray
        A 2D array (cells x sites). Each entry indicates the mutation state for a specific site in a specific cell.
    mut_prob : np.ndarray
        A 2D array (cells x sites). Each entry indicates the probability of a mutation at a site for a particular cell.
    beta : float
        A parameter for (1 + beta * t)^(-alpha) (scale factor in the integrals).
    max_iter : int, optional
        Maximum number of iterations for alternating optimization. Default is 100.
    tol : float, optional
        Convergence threshold for alpha_vector. Default is 1e-6.

    Returns
    -------
    alpha_vector : np.ndarray
        The optimized alpha_vector.
    MRCA_matrix : np.ndarray
        The optimized MRCA_matrix.
    """
    num_cells = sc_mat.shape[0]
    # Initialize alpha_vector (one alpha per site)
    alpha_vector = np.ones(sc_mat.shape[1])
    # Initialize MRCA_matrix as a zero matrix
    MRCA_matrix = np.zeros((num_cells, num_cells))

    for iteration in range(max_iter):
        print(f"Iteration {iteration + 1}...")

        # Step 1: Fix alpha_vector, optimize MRCA_matrix
        MRCA_matrix = compute_MRCA_time(t_S, sc_mat, mut_prob, alpha_vector, beta)

        # Step 2: Fix MRCA_matrix, optimize alpha_vector
        new_alpha_vector = compute_alpha(t_S, sc_mat, mut_prob, beta, MRCA_matrix)

        # Check convergence
        diff = np.linalg.norm(new_alpha_vector - alpha_vector)
        print(f"Difference in alpha_vector: {diff}")
        if diff < tol:
            print("Convergence achieved!")
            alpha_vector = new_alpha_vector
            break

        alpha_vector = new_alpha_vector

    return alpha_vector, MRCA_matrix



def tree_to_newick(
    node: ClusterNode,
    leaf_number: int,
    leaf_nodes: Dict[int, str],
    internal_nodes_map: Dict[int, str]
) -> str:
    """
    Convert a scipy ClusterNode to a partial Newick-formatted string (without the semicolon).

    Parameters
    ----------
    node : ClusterNode
        A node in the scipy hierarchical clustering tree.
    leaf_number : int
        The total number of leaf nodes (cells) in this tree.
    leaf_nodes : Dict[int, str]
        A dictionary mapping leaf node IDs to their names (e.g., {0: "Cell0", 1: "Cell1", ...}).
    internal_nodes_map : Dict[int, str]
        A dictionary mapping (node.id - leaf_number + 1) to internal node labels (e.g., {1: "Node1", 2: "Node2", ...}).

    Returns
    -------
    str
        A Newick-formatted subtree string for the given node.
    """
    # If the node is a leaf, return its corresponding label from leaf_nodes
    if node.is_leaf():
        return leaf_nodes[node.id]
    else:
        # Recursively build left and right subtrees
        left_subtree = tree_to_newick(node.left, leaf_number, leaf_nodes, internal_nodes_map)
        right_subtree = tree_to_newick(node.right, leaf_number, leaf_nodes, internal_nodes_map)

        # Retrieve the internal node label from internal_nodes_map
        # node.dist is the branch length from this internal node to its children
        return f"({left_subtree},{right_subtree}){internal_nodes_map[(node.id - leaf_number + 1)]}:{node.dist}"


def reconstrust_cell_lineage_tree(
    distance_matrix: np.ndarray,
    leaf_number: int,
    clustering_method: str = "average",
    leaf_nodes_dict: Dict[int, str] = None,
    internal_nodes_dict: Dict[int, str] = None
) -> Tree:
    """
    Build a cell lineage tree from a distance matrix using hierarchical clustering,
    convert it to a Newick string, then parse it into a dendropy Tree.

    Parameters
    ----------
    distance_matrix : np.ndarray
        A 2D numpy array representing pairwise distances (size: [leaf_number x leaf_number]).
    leaf_number : int
        The total number of leaf nodes/cells.
    clustering_method : str, optional
        The linkage method for SciPy's hierarchical clustering ("single", "complete", "average", "ward", etc.).
        Default is "average".
    leaf_nodes_dict : Dict[int, str], optional
        A mapping from leaf IDs to labels (e.g., {0: "Cell0", 1: "Cell1"}). If None, defaults to "Cell{i}".
    internal_nodes_dict : Dict[int, str], optional
        A mapping from (node.id - leaf_number + 1) to internal node labels.
        If None, defaults to {1: "Node1", 2: "Node2", ...}.

    Returns
    -------
    Tree
        A dendropy Tree object parsed from the generated Newick string.
    List[str]
        A list of internal node labels found in the tree.
    """
    # Perform hierarchical clustering
    cluster = linkage(distance_matrix, method=clustering_method)

    # Convert the linkage output to a ClusterNode tree
    t = to_tree(cluster, rd=False)

    # Default leaf node naming if not provided
    if leaf_nodes_dict is None:
        leaf_nodes_dict = {i: f"Cell{i}" for i in range(leaf_number)}

    # Default internal node naming if not provided
    if internal_nodes_dict is None:
        internal_nodes_dict = {i: f"IP{i}" for i in range(1, leaf_number)}

    # Convert the tree to Newick format (partial) and add the semicolon at the end
    newick_string = tree_to_newick(t, leaf_number, leaf_nodes_dict, internal_nodes_dict) + ";"

    # Parse the Newick string into a dendropy Tree
    cell_lineage_tree = dendropy.Tree.get(data=newick_string, schema='newick')

    # Collect internal node labels
    internal_nodes = []
    for node in cell_lineage_tree.postorder_node_iter():
        if node.label is not None:
            internal_nodes.append(node.label)

    return cell_lineage_tree, internal_nodes

def annotate_internal_nodes(tree: Tree, internal_nodes: List[str],cell_types: Dict[str,str]) -> Dict:
    """
    Annotate internal nodes in the given dendropy tree with the leaf labels that descend from each internal node.
    Also build an annotation dictionary mapping each internal node to the cell types of its descendant leaves.

    Parameters
    ----------
    tree : dendropy.Tree
        The dendropy tree object to annotate.
    internal_nodes : List[str]
        A list of internal node labels in the tree.

    Returns
    -------
    descendants : Dict
        A dictionary { internal_node_id: [list_of_leaf_labels] }.
    annotation : Dict
        A dictionary { internal_node_id: [list_of_cell_types] }.
    """
    descendants = {}
    # For each internal node label, find the corresponding node in the tree
    for internal_node_id in internal_nodes:
        for node in tree.postorder_node_iter():
            if node.label == internal_node_id and node.is_internal():
                internal_node = node
                break
        # Collect leaf nodes under this internal node
        subtree_leaf_nodes = [child for child in internal_node if child.is_leaf()]

        leaves = []
        for leaf in subtree_leaf_nodes:
            leaves.append(leaf.taxon.label)

        descendants.update({internal_node_id: leaves})

    annotation = descendants.copy()
    for key in descendants.keys():
        annotation[key] = []
        for value in descendants[key]:
            annotation[key].extend([cell_types[value]])
        annotation[key] = list(sorted(set(annotation[key])))

    return descendants, annotation


def find_child_nodes(tree: Tree, parent_node_label: str) -> List[str]:
    """
    Retrieve the labels of child nodes for a given parent node label within a dendropy Tree.
    If a child is a leaf (and thus has no 'label'), its taxon label is returned instead.

    Parameters
    ----------
    tree : dendropy.Tree
        A dendropy tree containing the parent node.
    parent_node_label : str
        Label of the parent node for which to find child nodes.

    Returns
    -------
    List[str]
        A list of child node labels or taxon labels (for leaves).
        If the parent node is not found or is a leaf itself, an empty list is returned.
    """
    parent_node = tree.find_node_with_label(parent_node_label)
    if parent_node and not parent_node.is_leaf():
        children = []
        for child in parent_node.child_nodes():
            if child.label is None:
                # If child is a leaf, use its taxon label
                children.append(child.taxon.label)
            else:
                children.append(child.label)
        return children
    return []


def average_separation_time(
    cell_type_1: str,
    cell_type_2: str,
    tree: Tree,
    cell_types: Dict[str,str],
    descendants: Dict[str,List[str]],
    internal_nodes: List[str]
):
    """
    Calculate the mean distance-from-root for internal nodes that distinctly separate
    two given cell types in the provided dendropy Tree.

    Parameters
    ----------
    cell_type_1 : str
        The first cell type of interest.
    cell_type_2 : str
        The second cell type of interest.
    tree : dendropy.Tree
        The dendropy tree to analyze.
    cell_types : dict
        A dictionary mapping leaf/taxon labels to their cell type.
    descendants : Dict
        A dictionary { internal_node_id: [list_of_leaf_labels] }
    internal_nodes : list of str
        Labels for the internal nodes in the tree.

    Returns
    -------
    float
        The average distance-from-root for all internal nodes that split cell_type_1 and cell_type_2
        into different child branches. If no such node exists, returns NaN.
    """
    separation_nodes = []
    separation_time = []

    # Check each internal node to see if it separates cell_type_1 and cell_type_2
    for internal_node_id in internal_nodes:
        children = find_child_nodes(tree, internal_node_id)
        descendant_celltypes = dict({children[0]: [], children[1]: []})

        for child in children:
            if child not in internal_nodes:
                if len([cell_types[child]]) == 1:
                    descendant_celltypes[child].append(cell_types[child])
                else:
                    descendant_celltypes[child].extend(cell_types[child])
            else:
                leaves = descendants[child]
                for leaf in leaves:
                    descendant_celltypes[child].append(cell_types[leaf])

        # Check whether this internal node indeed separates the two cell types distinctly
        if (
            (
                (cell_type_1 in list(descendant_celltypes.values())[0])
                or (cell_type_1 in list(descendant_celltypes.values())[1])
            )
            and
            (
                (cell_type_2 in list(descendant_celltypes.values())[0])
                or (cell_type_2 in list(descendant_celltypes.values())[1])
            )
        ) and (
            (
                (cell_type_1 in list(descendant_celltypes.values())[0])
                ^ (cell_type_2 in list(descendant_celltypes.values())[0])
            )
            and
            (
                (cell_type_1 in list(descendant_celltypes.values())[1])
                ^ (cell_type_2 in list(descendant_celltypes.values())[1])
            )
        ):
            separation_nodes.append(internal_node_id)

    # For each node that separates the two cell types, compute distance_from_root()
    for separation_node_id in separation_nodes:
        for node in tree.postorder_node_iter():
            if node.label == separation_node_id:
                target_node = node
                break
        separation_time.append(target_node.distance_from_root())

    return np.mean(separation_time)


def average_separation_time_dict(tree: Tree, cell_types: Dict[int, str], descendants: Dict[str,List[str]], internal_nodes: List[str]) -> Dict[tuple, float]:
    """
    Build a dictionary of pairwise separation times for every unique pair of cell types.

    Parameters
    ----------
    tree : dendropy.Tree
        The dendropy tree containing the nodes.
    cell_types : dict
        A dictionary mapping leaf/taxon labels to their cell types.
    descendants : Dict
        A dictionary { internal_node_id: [list_of_leaf_labels] }
    internal_nodes : list of str
        Labels for the internal nodes in the tree.
    
    Returns
    -------
    Dict[tuple, float]
        A dictionary where each key is a tuple (cell_type_1, cell_type_2) and
        each value is the average separation time between these two types,
        sorted in descending order by time.
    """
    all_types = list(set(cell_types.values()))
    time_dict = {}
    for pair in combinations(all_types, 2):
        time = average_separation_time(pair[0], pair[1], tree, cell_types, descendants, internal_nodes)
        time_dict.update({pair: time})

    return dict(sorted(time_dict.items(), key=lambda item: item[1], reverse=True))


class UPGMANode:
    """
    A Node class used for building the UPGMA tree.
    Attributes:
        name (str): Name of the node (leaf name or internal node name like "Node1").
        children (list): List of child nodes (empty if leaf).
        height (float): The distance from this node down to the leaves in the UPGMA context.
    """
    def __init__(self, name, children=None, height=0.0):
        self.name = name
        self.children = children if children else []
        self.height = height

    def is_leaf(self):
        """Check if the node is a leaf."""
        return len(self.children) == 0


def upgma(distance_df):
    """
    Perform UPGMA on the given distance matrix (pandas DataFrame).
    Rows and columns must represent the same set of leaf labels.
    
    Args:
        distance_df (pd.DataFrame): Symmetric distance matrix for leaves.
        
    Returns:
        UPGMANode: The root node of the constructed UPGMA tree.
    """
    # Each initial cluster is a single leaf node
    clusters = {
        name: (UPGMANode(name=name), 1)    # (node_object, cluster_size)
        for name in distance_df.index
    }

    # Copy the distance matrix to update dynamically
    cur_df = distance_df.copy()

    internal_node_count = 0

    # While more than one cluster remains, merge the closest pair
    while len(clusters) > 1:
        # 1) Find the two clusters with the smallest distance
        min_dist = float('inf')
        pair_to_merge = (None, None)

        names_list = list(clusters.keys())
        for i in range(len(names_list)):
            for j in range(i+1, len(names_list)):
                n1, n2 = names_list[i], names_list[j]
                d = cur_df.loc[n1, n2]
                if d < min_dist:
                    min_dist = d
                    pair_to_merge = (n1, n2)

        # 2) Merge these two clusters
        n1, n2 = pair_to_merge
        (node1, size1) = clusters[n1]
        (node2, size2) = clusters[n2]

        internal_node_count += 1
        new_node_name = f"IP{internal_node_count}"

        # In UPGMA, the new internal node's height is half the distance
        new_height = min_dist / 2.0
        new_node = UPGMANode(
            name=new_node_name,
            children=[node1, node2],
            height=new_height
        )

        # Combined cluster size
        new_size = size1 + size2

        # Remove merged clusters, add new cluster
        del clusters[n1]
        del clusters[n2]
        clusters[new_node_name] = (new_node, new_size)

        # 3) Update the distance matrix for the newly formed cluster
        for other_name in clusters.keys():
            if other_name == new_node_name:
                continue
            (other_node, other_size) = clusters[other_name]

            # UPGMA average distance formula
            dist_new_other = (
                size1 * cur_df.loc[n1, other_name] +
                size2 * cur_df.loc[n2, other_name]
            ) / new_size

            cur_df.loc[new_node_name, other_name] = dist_new_other
            cur_df.loc[other_name, new_node_name] = dist_new_other

        # Remove rows and columns of the merged clusters
        cur_df = cur_df.drop(index=[n1, n2], columns=[n1, n2])

    # Only one cluster (the root) remains
    root_name = list(clusters.keys())[0]
    root_node, _ = clusters[root_name]
    return root_node



def to_newick(node):
    """
    Convert a UPGMA tree node (root) to Newick string.
    The length after the colon is the local 'height' used for visualization.
    
    Args:
        node (UPGMANode): The root of the UPGMA tree.
    
    Returns:
        str: Newick-format string representing the tree.
    """
    if node.is_leaf():
        return f"{node.name}:{node.height:.2f}"
    else:
        children_str = ",".join([to_newick(child) for child in node.children])
        return f"({children_str}){node.name}:{node.height:.2f}"


class ParsedNode:
    """
    A Node class used for the parsed Newick tree, where `height` is interpreted
    as the distance from this node down to the leaves (common in UPGMA trees).
    """
    def __init__(self, name='', height=0.0):
        self.name = name
        self.height = height
        self.children = []

    def __repr__(self):
        return f"ParsedNode(name={self.name}, height={self.height}, children={len(self.children)})"



class NewickParser:
    """
    Simple Newick parser: 
    - Assumes each ':number' represents 'distance from this node to the leaves'.
    - Ignores any bootstrap/comment notations.
    - Expected format: (child1,child2)NodeName:Height  or LeafName:Height
    """
    def __init__(self, newick_str):
        self.s = newick_str.strip()
        if self.s.endswith(';'):
            self.s = self.s[:-1]
        self.pos = 0

    def current_char(self):
        if self.pos >= len(self.s):
            return ''
        return self.s[self.pos]

    def advance(self):
        self.pos += 1

    def parse_subtree(self):
        """
        Parse a subtree which might be:
          (child1,child2)NodeX:Height
          or
          LeafName:Height
        """
        node = ParsedNode()

        # If '(' => internal node with children
        if self.current_char() == '(':
            self.advance()
            node.children = self.parse_children()
            if self.current_char() == ')':
                self.advance()
            else:
                raise ValueError("Missing closing parenthesis")

        # Then node name
        node.name = self.parse_node_name()

        # Then optional :Height
        if self.current_char() == ':':
            self.advance()
            node.height = self.parse_float()

        return node

    def parse_children(self):
        """
        Parse comma-separated list of subtrees (children)
        """
        children = []
        while True:
            child = self.parse_subtree()
            children.append(child)

            if self.current_char() == ',':
                self.advance()
                continue
            else:
                break
        return children

    def parse_node_name(self):
        """
        Parse node name (alphanumeric/underscore)
        """
        pattern = re.compile(r'[A-Za-z0-9_]+')
        match = pattern.match(self.s, self.pos)
        if match:
            name = match.group(0)
            self.pos = match.end()
            return name
        return ''

    def parse_float(self):
        """
        Parse a floating-point number (possibly scientific notation)
        """
        pattern = re.compile(r'[+-]?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?')
        match = pattern.match(self.s, self.pos)
        if not match:
            raise ValueError("Expected a floating-point number")
        val = float(match.group(0))
        self.pos = match.end()
        return val


def parse_newick_to_tree(newick_str):
    """
    Parse a UPGMA-styled Newick string into a tree of ParsedNode objects.
    
    Args:
        newick_str (str): Newick-formatted string (should end in ';').
    
    Returns:
        ParsedNode: The root node of the parsed tree.
    """
    parser = NewickParser(newick_str)
    root = parser.parse_subtree()

    if parser.pos < len(parser.s):
        raise ValueError("Extra unparsed characters at the end of Newick string.")
    return root


def build_edge_dict_upgma_normalized(root):
    """
    Given a parsed tree (root), where each node's 'height' is 
    the distance from this node down to the leaves, compute the 
    normalized edge length for every Parent->Child.
    
    For each Parent->Child:
      edge_length = (parent.height - child.height) / root.height
    
    Args:
        root (ParsedNode): Root of the parsed Newick tree.
    
    Returns:
        dict: Key: "ParentName->ChildName", Value: normalized edge length (float).
    """
    edge_dict = {}
    total_height = root.height  # The distance from root to leaves

    def dfs(parent, node):
        if parent is not None:
            raw_length = parent.height - node.height
            # Normalize by the total height at the root
            if total_height > 0:
                norm_length = raw_length / total_height
            else:
                norm_length = 0.0

            key = f"{parent.name}->{node.name}"
            edge_dict[key] = norm_length

        for ch in node.children:
            dfs(node, ch)

    dfs(None, root)
    return edge_dict


def upgma_to_edge_dict(distance_df):
    """
    High-level function that:
      1. Builds a UPGMA tree from the input distance matrix.
      2. Converts the tree to Newick.
      3. Parses that Newick to interpret node heights as distance-to-leaves.
      4. Builds a dictionary of normalized edge lengths (0.0 to 1.0).
    
    Args:
        distance_df (pd.DataFrame): Input distance matrix (symmetric).
    
    Returns:
        dict: A dictionary where each key is "ParentName->ChildName" and
              each value is the normalized edge length between them.
    """
    root_upgma = upgma(distance_df)
    newick_str = to_newick(root_upgma) + ";"
    parsed_root = parse_newick_to_tree(newick_str)
    edge_dict = build_edge_dict_upgma_normalized(parsed_root)

    return edge_dict

def parse_edge_dict(edge_dict:Dict[str,float]):
    edges = []
    for key in edge_dict:
        start, end = key.split('->')
        weight = edge_dict[key]
        fate_map_edge = FateMapEdge(start,end,weight)
        edges.append(fate_map_edge)
    return edges

def cell_types_distance(MRCA_matrix, cell_types):
    """
    Construct a cell-type distance matrix from an MRCA (Most Recent Common Ancestor) matrix.

    This function:
      1. Interprets the MRCA_matrix as a pairwise distance among cells.
      2. Reconstructs a cell lineage tree using hierarchical clustering.
      3. Annotates the internal nodes of the tree with their descendant leaves.
      4. Calculates the average separation time for each pair of cell types.
      5. Builds a final NxN distance matrix (where N is the number of unique cell types)
         and returns it as a pandas DataFrame.

    Parameters
    ----------
    MRCA_matrix : np.ndarray
        A square matrix representing pairwise MRCA distances among cells.
        Its shape is [num_cells x num_cells].
    cell_types : dict
        A dictionary mapping leaf (cell) labels to their corresponding cell types.
        For example: {"Cell0": "TypeA", "Cell1": "TypeB", ...}.

    Returns
    -------
    distance_df : pd.DataFrame
        A symmetric DataFrame where both rows and columns represent unique cell types.
        Each entry is the average separation distance (or time) between two cell types.
    """
    # Determine the number of leaf nodes based on the MRCA_matrix dimensions
    leaf_number = len(MRCA_matrix)

    # Reconstruct a cell lineage tree using hierarchical clustering (average linkage by default)
    cell_lineage_tree, internal_nodes = reconstrust_cell_lineage_tree(
        MRCA_matrix, 
        leaf_number, 
        clustering_method="average"
    )

    # Annotate the internal nodes of the tree to figure out which leaves belong to them
    descendants, annotation = annotate_internal_nodes(cell_lineage_tree, internal_nodes,cell_types)

    # Calculate the average separation time for each pair of cell types
    separation_time = average_separation_time_dict(
        cell_lineage_tree,
        cell_types,
        descendants,
        internal_nodes
    )

    # Build a distance matrix among the unique cell types
    nodes = set()
    for key in separation_time.keys():
        # Each key is a tuple of (cell_type1, cell_type2)
        nodes.update(key)

    # Sort the list of cell types for consistent ordering
    nodes = sorted(nodes)
    size = len(nodes)

    # Initialize a square matrix of shape (size x size)
    distance_matrix = np.zeros((size, size))

    # Create a dictionary to map each cell type to its row/column index
    node_index = {node: idx for idx, node in enumerate(nodes)}

    # Populate the matrix with separation distances
    for (node1, node2), distance in separation_time.items():
        idx1 = node_index[node1]
        idx2 = node_index[node2]
        distance_matrix[idx1, idx2] = distance
        distance_matrix[idx2, idx1] = distance

    # Convert to a pandas DataFrame for easy viewing/manipulation
    distance_df = pd.DataFrame(distance_matrix, index=nodes, columns=nodes)

    return distance_df


def upgma_to_edge_dict(distance_df):
    """
    Builds a UPGMA tree from the input distance matrix, then generates a dictionary
    of normalized edge lengths in the format { "Parent->Child": normalized_length }.

    Steps:
      1. Perform UPGMA on distance_df to get a tree.
      2. Convert the resulting tree to a Newick string.
      3. Parse the Newick string to interpret node heights as distance-to-leaves.
      4. Build a dictionary of normalized edge lengths, where each edge length is
         (parent.height - child.height) / root.height.

    Parameters
    ----------
    distance_df : pd.DataFrame
        A symmetric distance matrix where rows and columns represent entities
        to be clustered.

    Returns
    -------
    dict
        A dictionary where each key is "ParentName->ChildName" and
        each value is the normalized edge length (float) between them.
    """
    # Step 1: Build a UPGMA tree
    root_upgma = upgma(distance_df)

    # Step 2: Convert that tree to a Newick string
    newick_str = to_newick(root_upgma) + ";"

    # Step 3: Parse the Newick string (treat node height as distance to leaves)
    parsed_root = parse_newick_to_tree(newick_str)

    # Step 4: Build a dictionary of normalized edge lengths
    edge_dict = build_edge_dict_upgma_normalized(parsed_root)
    return edge_dict

def load_fate_map_topology(fate_map_path):
    """
    Load a dictionary from a CSV file.

    Parameters:
        fate_map_path (str): Path to the input CSV file.

    Returns:
        dict: The reconstructed dictionary.
    """
    df = pd.read_csv(fate_map_path)
    edge_dict = {f"{row['Parent']}->{row['Child']}": row['Length'] for _, row in df.iterrows()}
    return edge_dict

def parse_edge_dict(edge_dict: Dict[str, float]):
    """
    Convert an edge dictionary (key = "Parent->Child", value = weight)
    into a list of FateMapEdge objects.

    Parameters
    ----------
    edge_dict : Dict[str, float]
        A dictionary where each key has the format "ParentName->ChildName"
        and the value is the edge weight or distance.

    Returns
    -------
    list
        A list of FateMapEdge objects constructed from the key-value pairs.
    """
    edges = []
    for key, weight in edge_dict.items():
        start, end = key.split('->')
        grn_edge = FateMapEdge(start, end, weight)
        edges.append(grn_edge)
    return edges


def construct_fate_map(sc_mat, cell_types, t_S, beta=1, max_iter=100, tol=1e-2, output_edges=False):
    """
    Construct a fate map based on single-cell matrix ``sc_mat``, cell types, 
    and inferred alpha/MRCA matrices. This is a high-level pipeline.

    **Steps:**

    1. **Compute default allele probabilities.**
    2. **Perform alternating optimization to get final alpha and MRCA_matrix.**
    3. **Build a cell types distance DataFrame using ``cell_types_distance``.**
    4. **Convert that distance DataFrame into a dictionary of normalized UPGMA edges.**
    5. **Parse the edge dictionary into ``GRNEdge`` objects.**
    6. **Wrap edges into a ``FateMap`` object.**

    :param np.ndarray sc_mat: 
        Single-cell matrix, where rows typically represent cells 
        and columns represent sites.
    :param dict cell_types: 
        Mapping from cell labels to their cell type names.
    :param float t_S: 
        A sampling time or similar parameter for the optimization.
    :param float beta: 
        A parameter for the (1 + beta * t) model in the integrals. 
        **Default:** 1.
    :param int max_iter: 
        Maximum number of iterations for alternating optimization. 
        **Default:** 1000.
    :param float tol: 
        Convergence threshold. 
        **Default:** 1e-5.
    :param bool output_edges: 
        If True, returns a tuple of (edge, fate_map). 
        Otherwise returns just the fate_map.

    :returns: 
        - If ``output_edges`` is False: 
          
          A ``FateMap`` object containing the inferred edges 
          for further analysis/visualization.
        
        - If ``output_edges`` is True: 
          
          A tuple ``(edge, fate_map)`` where ``edge`` represents 
          the GRNEdge objects.
    :rtype: FateMap or tuple
    """

    # 1. Compute default allele probabilities
    mut_prob = compute_default_allele_prob_new(sc_mat)

    # 2. Perform alternating optimization to get alpha and MRCA
    final_alpha, MRCA_matrix = alternating_optimization(t_S, sc_mat, mut_prob, beta, max_iter, tol)

    # 3. Build a cell-type distance DataFrame
    distance_df = cell_types_distance(MRCA_matrix, cell_types)

    # 4. Convert distance DataFrame to edge dict via UPGMA
    edge_dict = upgma_to_edge_dict(distance_df)

    # 5. Parse edge dict to create GRNEdge objects
    edge = parse_edge_dict(edge_dict)

    # 6. Construct the FateMap (wrap the edges)
    fate_map = FateMap(edge)
    
    if output_edges:
        return edge, fate_map
    else:
        return fate_map

