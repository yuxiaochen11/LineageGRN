def get_regulators_for_target_gene(target_gene_id: str, node_id: str, input_path: str, regulator_names: list[str]):
    """
    Retrieve regulator information for a specific target gene.

    :param target_gene_id: **The ID of the target gene.**
    :type target_gene_id: str
    :param node_id: **The ID of the node.**
    :type node_id: str
    :param input_path: **The path to the input files.**
    :type input_path: str
    :param regulator_names: **A list of names of regulators.**
    :type regulator_names: list of str

    :returns: **A dictionary** containing the regulator information, with keys as regulator names
              and values as their corresponding lambda values.
    :rtype: dict
    """
    file_path = input_path + '/' + 'target_gene_' + target_gene_id + '.csv'
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

    :param regulator_id: **The ID of the regulator.**
    :type regulator_id: str
    :param node_id: **The ID of the node.**
    :type node_id: str
    :param input_path: **The path to the input files.**
    :type input_path: str
    :param regulator_names: **A list of names of regulators.**
    :type regulator_names: list of str

    :returns: **A dictionary** where the keys are target gene IDs and the values are
              the corresponding lambda values for the specified regulator.
    :rtype: dict
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



def identify_key_genes_differentiation(input_path, fate_map, threshold, regulator_names, target_gene_names, ancestor_node_id):
    """
    Identify key regulators involved in cell fate decision events.

    :param input_path: **The path to the input files** containing regulatory data.
    :type input_path: str
    :param fate_map: **An object** that contains information about the nodes and edges.
    :type fate_map: object
    :param threshold: **A threshold value** to filter regulatory effects.
    :type threshold: float
    :param regulator_names: **A list of names** of regulators.
    :type regulator_names: list of str
    :param target_gene_names: **A list of names** of target genes.
    :type target_gene_names: list of str
    :param ancestor_node_id: **The ID of the ancestor node** from which to analyze child nodes.
    :type ancestor_node_id: str

    :returns: **A DataFrame** containing key regulators and their associated metrics for the given ancestor node.
    :rtype: pd.DataFrame
    """
    # Get child nodes from the ancestor node
    child_nodes = [fate_map.nodes[ancestor_node_id].directed_edges[i].end for i in range(2)]
    
    # Combine high expression target genes from the child nodes
    high_expression_target_genes_in_children = (
        fate_map.nodes[child_nodes[0]].high_expression_genes_in_leaf + 
        fate_map.nodes[child_nodes[1]].high_expression_genes_in_leaf
    )
    
    # Identify high expression target genes present in both target genes and child nodes
    high_expression_target_genes_in_child_nodes = list(set(target_gene_names) & set(high_expression_target_genes_in_children))
    
    # Retrieve the GRN for the ancestor node
    grn_dict = get_dynamic_networks(input_path, fate_map, threshold, regulator_names, target_gene_names)

    # Create a DataFrame for regulators
    regulators_df = pd.DataFrame(regulator_names, columns=['regulator_id'])
    regulators_df.index = regulators_df['regulator_id']

    # Calculate the sum of regulatory values for high expression target genes
    column_sums = grn_dict[ancestor_node_id].loc[high_expression_target_genes_in_child_nodes, :].sum(axis=0)
    grn_df = grn_dict[ancestor_node_id].loc[high_expression_target_genes_in_child_nodes, :][column_sums[column_sums > 0].index]

    # Count the number of positive and negative regulatory effects
    positive_regulation_number = pd.DataFrame(grn_df.apply(lambda x: (x > 0).sum()))
    negative_regulation_number = pd.DataFrame(grn_df.apply(lambda x: (x < 0).sum()))

    # Sum regulatory values across all target genes
    regulatory_value_sum = pd.DataFrame(grn_df.sum(axis=0))
    
    # Merge dataframes to create a summary of key regulators
    key_regulators_df = pd.merge(regulatory_value_sum, positive_regulation_number, left_index=True, right_index=True, how='outer')
    key_regulators_df = pd.merge(key_regulators_df, negative_regulation_number, left_index=True, right_index=True, how='outer')

    # Rename columns for clarity
    key_regulators_df.columns = ['positive_regulatory_strength', 'positive_regulation_number', 'negative_regulation_number']
    
    # Calculate probability and activation scale
    key_regulators_df['DRS'] = key_regulators_df['positive_regulatory_strength'].div(key_regulators_df['positive_regulation_number'])
    key_regulators_df['DRC'] = key_regulators_df['positive_regulation_number'].div(
        key_regulators_df['positive_regulation_number'] + key_regulators_df['negative_regulation_number']
    )
    key_regulators_df = key_regulators_df.drop_duplicates()  # Remove duplicates
    
    # Join with regulators DataFrame and sort by regulatory metrics
    key_regulators_df = regulators_df.join(key_regulators_df, how='left')
    key_regulators_df = key_regulators_df.sort_values(
        by=['negative_regulation_number', 'positive_regulation_number', 'DRC'], 
        ascending=[True, False, False]
    )
    
    key_regulators_df['node_id'] = [ancestor_node_id] * key_regulators_df.shape[0]  # Add ancestor node ID to the DataFrame

    return key_regulators_df


def identify_key_genes_fate_bias(grn_dict, lineage, Tar_1, Tar_2, regulator_names):
    """
    Analyze fate bias genes based on the regulatory network data.

    :param grn_dict: **A dictionary** containing the gene regulatory networks (GRNs) for different lineages.
    :type grn_dict: dict
    :param lineage: **The lineage** to be analyzed.
    :type lineage: str
    :param Tar_1: **The name of the first target gene**.
    :type Tar_1: str
    :param Tar_2: **The name of the second target gene**.
    :type Tar_2: str
    :param regulator_names: **A list of names** of regulators to analyze.
    :type regulator_names: list of str

    :returns: **Two DataFrames** containing key regulators for the specified target genes, sorted by regulatory bias.
    :rtype: tuple of pd.DataFrame
    """
    regulation_data = []
    grn_df = grn_dict[lineage]
    
    for regulator in regulator_names:
        regulator_values = grn_df[regulator]

        tar_1_values = regulator_values[Tar_1]
        tar_2_values = regulator_values[Tar_2]

        positive_tar_1_count = sum(tar_1_values > 0)
        negative_tar_2_count = sum(tar_2_values < 0)
        positive_tar_1_strength = tar_1_values[tar_1_values > 0].mean() if positive_tar_1_count > 0 else 0
        positive_tar_2_count = sum(tar_2_values > 0)
        negative_tar_1_count = sum(tar_1_values < 0)
        positive_tar_2_strength = tar_2_values[tar_2_values > 0].mean() if positive_tar_2_count > 0 else 0
        
        regulation_data.append({
            'Regulatory gene': regulator,
            'Tar_1_positive_count': positive_tar_1_count,
            'Tar_2_negative_count': negative_tar_2_count,
            'Tar_1_strength': positive_tar_1_strength,
            'Tar_2_positive_count': positive_tar_2_count,
            'Tar_1_negative_count': negative_tar_1_count,
            'Tar_2_strength': positive_tar_2_strength
        })
    
    key_regulators_df = pd.DataFrame(regulation_data)

    # Calculate regulatory bias for both target genes
    key_regulators_of_T1 = key_regulators_df.copy()
    key_regulators_of_T2 = key_regulators_df.copy()

    # Calculate Positive Regulatory Bias, Negative Regulatory Bias, and Strength Bias for Tar_1 and Tar_2
    key_regulators_of_T1['FAP'] = (
        key_regulators_of_T1['Tar_1_positive_count'] 
        / (key_regulators_of_T1['Tar_1_positive_count'] + key_regulators_of_T1['Tar_2_positive_count'])
    )
    key_regulators_of_T1['FRP'] = (
        key_regulators_of_T1['Tar_2_negative_count'] 
        / (key_regulators_of_T1['Tar_2_negative_count'] + key_regulators_of_T1['Tar_1_negative_count'])
    )
    key_regulators_of_T1['FBS'] = (
        key_regulators_of_T1['Tar_1_strength'] 
        / (key_regulators_of_T1['Tar_1_strength'] + key_regulators_of_T1['Tar_2_strength'])
    )
    key_regulators_of_T1['FBP'] = key_regulators_of_T1['FAP'] + key_regulators_of_T1['FRP']
    key_regulators_of_T1 = key_regulators_of_T1.fillna(0)

    key_regulators_of_T2['FAP'] = (
        key_regulators_of_T2['Tar_2_positive_count'] 
        / (key_regulators_of_T2['Tar_2_positive_count'] + key_regulators_of_T2['Tar_1_positive_count'])
    )
    key_regulators_of_T2['FRP'] = (
        key_regulators_of_T2['Tar_1_negative_count'] 
        / (key_regulators_of_T2['Tar_1_negative_count'] + key_regulators_of_T2['Tar_2_negative_count'])
    )
    key_regulators_of_T2['FRP'] = key_regulators_of_T2['FRP'].fillna(0)
    key_regulators_of_T2['FBS'] = (
        key_regulators_of_T2['Tar_2_strength'] 
        / (key_regulators_of_T2['Tar_2_strength'] + key_regulators_of_T2['Tar_1_strength'])
    )
    key_regulators_of_T2['FBP'] = key_regulators_of_T2['FAP'] + key_regulators_of_T2['FRP']
    key_regulators_of_T2 = key_regulators_of_T2.fillna(0)

    # Split the data into subsets based on regulatory strength
    subsetA_T1 = key_regulators_of_T1[key_regulators_of_T1['FBS'] <= 0.5].copy()
    subsetB_T1 = key_regulators_of_T1[key_regulators_of_T1['FBS'] > 0.5].copy()

    subsetA_T1 = subsetA_T1.sort_values(by='FBP', ascending=True)
    subsetB_T1 = subsetB_T1.sort_values(by=['FBP','FBS','FRP','FAP'], 
                                        ascending=[True,False,True,True])
    key_regulators_of_T1 = pd.concat([subsetA_T1, subsetB_T1], ignore_index=True)
    
    subsetA_T2 = key_regulators_of_T2[key_regulators_of_T2['FBS'] > 0.5].copy()
    subsetB_T2 = key_regulators_of_T2[key_regulators_of_T2['FBS'] <= 0.5].copy()

    key_regulators_of_T2 = pd.concat([subsetA_T2, subsetB_T2], ignore_index=True)

    return key_regulators_of_T1, key_regulators_of_T2


def cluster_regulatory_interactions(input_path, fate_map, threshold, regulator_names, target_gene_names, regulator_number, target_number, n_clusters, m, max_iter=100, theta=1e-5, seed=0):
    """
    Run the Fuzzy C-Means clustering algorithm.

    :param input_path: **The path to the input files** containing regulatory data.
    :type input_path: str
    :param fate_map: **An object** that contains information about nodes and edges.
    :type fate_map: object
    :param threshold: **A threshold** to filter regulatory interactions.
    :type threshold: float
    :param regulator_names: **A list of names** of regulators.
    :type regulator_names: list of str
    :param target_gene_names: **A list of names** of target genes.
    :type target_gene_names: list of str
    :param regulator_number: **The number of regulators**.
    :type regulator_number: int
    :param target_number: **The number of target genes**.
    :type target_number: int
    :param n_clusters: **The number of clusters** to form.
    :type n_clusters: int
    :param m: **Fuzziness parameter** (typically greater than 1).
    :type m: float
    :param max_iter: **Maximum number of iterations** to run the algorithm (default is 100).
    :type max_iter: int, optional
    :param theta: **Convergence threshold** for the weight change (default is 1e-5).
    :type theta: float, optional
    :param seed: **Random seed** for reproducibility (default is 0).
    :type seed: int, optional

    :returns: **A tuple** containing:
        - DataFrame: The input data matrix with clustering information.
        - DataFrame: The coordinates of the cluster centers.
        - DataFrame: The membership weight matrix.
    :rtype: tuple of pd.DataFrame
    """
    rng = np.random.RandomState(seed)  # Set the random seed
    grns_dict = get_dynamic_networks(input_path, fate_map, threshold, regulator_names, target_gene_names)
    
    # Initialize data matrix X
    X = pd.DataFrame(0, columns=list(grns_dict.keys()), index=list(range(regulator_number * target_number)))
    
    # Populate the data matrix based on regulatory interactions
    for node in list(grns_dict.keys()):
        grn = grns_dict[node]
        for regulator in range(regulator_number):
            for target in range(target_number):
                i = regulator * target_number + target
                if abs(grn.iloc[target, regulator]) > threshold:
                    X.loc[i, node] = 1
                else:
                    X.loc[i, node] = 0
    
    N, D = np.shape(X)
    weight = rng.uniform(size=(N, n_clusters))  # Initialize weights randomly
    weight = weight / np.sum(weight, axis=1, keepdims=True)  # Normalize weights

    # Main iteration loop for updating weights and centers
    for i in range(max_iter):
        weight_old = weight.copy()  # Save old weights for convergence check
        centers = get_centers(weight, X, m)  # Update centers
        weight = get_weight(X, centers, m)  # Update weights
        if np.linalg.norm(weight - weight_old) < theta:  # Check for convergence
            break
    
    # Create edge ID mapping for results
    edge_id_dict = {}
    for i in range(target_number * regulator_number):
        edge_id_dict.update({i: regulator_names[i // target_number] + '->' + target_gene_names[i % target_number]})
    
    X = pd.DataFrame(X)
    X.columns = list(grns_dict.keys())
    centers = pd.DataFrame(centers)
    weight_matrix = pd.DataFrame(weight)
    weight_matrix.columns = ['EdgeCluster_' + str(i) for i in range(1, (n_clusters + 1))]
    weight_matrix.index = list(edge_id_dict.values())
    
    return X, centers, weight_matrix


def identify_regulatory_interactions_specificity(input_path, fate_map, threshold, regulator_names, target_gene_names, n_clusters, weight_threshold, X, regulator_number, target_number, weight_matrix):
    """
    Map edge_clusters to nodes based on the regulatory interactions and weights.

    :param input_path: **The path to the input files** containing regulatory data.
    :type input_path: str
    :param fate_map: **An object** that contains information about nodes and edges.
    :type fate_map: object
    :param threshold: **A threshold** to filter regulatory interactions.
    :type threshold: float
    :param regulator_names: **A list of names** of regulators.
    :type regulator_names: list of str
    :param target_gene_names: **A list of names** of target genes.
    :type target_gene_names: list of str
    :param n_clusters: **The number of clusters** formed during clustering.
    :type n_clusters: int
    :param weight_threshold: **A threshold for considering the significance of weights**.
    :type weight_threshold: float
    :param X: **A DataFrame** representing the regulatory interactions.
    :type X: pd.DataFrame
    :param regulator_number: **The number of regulators**.
    :type regulator_number: int
    :param target_number: **The number of target genes**.
    :type target_number: int
    :param weight_matrix: **A DataFrame** containing membership weights for each edge and cluster.
    :type weight_matrix: pd.DataFrame

    :returns: **A DataFrame** mapping each cluster to nodes, showing the presence of edge_clusters in nodes.
    :rtype: pd.DataFrame
    """
    # Retrieve the gene regulatory networks for each node
    grns_dict = get_dynamic_networks(input_path, fate_map, threshold, regulator_names, target_gene_names)
    
    # Initialize a DataFrame to store the count of edge_clusters in each node
    edge_cluster_to_nodes = pd.DataFrame(0, columns=list(grns_dict.keys()), index=list(range(n_clusters)))
    
    # Count the number of edge_clusters associated with each node for each cluster
    for cluster in range(n_clusters):
        for edge in range(target_number * regulator_number):
            node_list = list(X.columns[(X.iloc[edge] == 1)])  # Get nodes connected by the edge
            if weight_matrix.iloc[edge, cluster] > weight_threshold:  # Check if weight exceeds threshold
                for j in node_list:
                    edge_cluster_to_nodes.loc[cluster, j] += 1  # Increment count for the node

    # Create a binary copy of the weight matrix based on the weight threshold
    weight_copy = weight_matrix.copy()
    weight_copy[weight_copy >= weight_threshold] = 1
    weight_copy[weight_copy < weight_threshold] = 0
    
    # Create a DataFrame to locate edge_clusters in nodes
    locate_edge_cluster_to_nodes_df = pd.DataFrame()
    for cluster in range(n_clusters):
        # Normalize the counts by the sum of weights for the respective cluster
        row = edge_cluster_to_nodes.iloc[cluster, :] / list(weight_copy.sum())[cluster]
        locate_edge_cluster_to_nodes_df = pd.concat([locate_edge_cluster_to_nodes_df, row], axis=1)
    
    # Set the column names to the weight DataFrame's columns
    locate_edge_cluster_to_nodes_df.columns = weight_matrix.columns
    
    return locate_edge_cluster_to_nodes_df



