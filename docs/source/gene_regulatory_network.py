class GRNInference:
    def __init__(self, atac_file_path, expression_file_path, fate_map, saved_dir):
        self.fate_map = fate_map
        self.atac_data = ATACData(atac_file_path)
        self.rna_data = RNAData(expression_file_path)
        self.saved_dir = saved_dir

    def _likelihood(self, leaves_grn_data, regulator_data, target_data):
        G = np.reshape(leaves_grn_data, newshape=[len(self.fate_map.node_leaves), len(self.rna_data.regulator_genes)])
        G = G[..., np.newaxis]
        D = target_data - np.einsum('ijk,ikl -> ijl', regulator_data, G)
        exp_sum = np.sum(-0.5 * SIGMA * np.einsum('ijk,ijk -> i', D, D)) / len(self.rna_data.cell_ids)
        constants_sum = -0.5 * len(self.rna_data.regulator_genes) * np.log(2 * np.pi) * len(self.fate_map.node_leaves) 
        result = exp_sum + constants_sum
        return result * LIKELIHOOD_WEIGHT

    def _atac_prior(self, leaves_grn_data, mdata):
        leaves_grn_data = np.stack((abs(leaves_grn_data) > 0.1, abs(leaves_grn_data) <= 0.1), axis=-1, dtype=np.int32)
        G = np.reshape(leaves_grn_data, newshape=[len(self.fate_map.node_leaves), len(self.rna_data.regulator_genes), 2])
        p = np.sum(np.einsum('ijk,ijk -> ij', mdata, G))
        return p * ATAC_PRIOR_WEIGHT

    def _lineage_prior(self, leaves_grn_data: 'np.ndarray', lambda_1: 'float', root_grn_data: 'np.ndarray', ancestor_nodes_space: 'np.ndarray', weights: 'np.ndarray'):
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
    
        grn_values = {}
        x = list(self.infer_leaf_grns(target_gene_id))

        leaf_node_ids = self.fate_map.node_leaves
        for idx, node_id in enumerate(leaf_node_ids):
            grn_value = x[1+idx*len(self.rna_data.regulator_genes):1+(idx+1)*len(self.rna_data.regulator_genes)]
            grn_values[node_id] = x[:1] + grn_value

        grn_values = self.infer_ancestor_grns(grn_values, target_gene_id)
        if saved:
            self._save(grn_values, target_gene_id)
            logger.info(f"Saved grn values for target_gene_id:{target_gene_id}")

        return grn_values

    def _reshape_likehood_input(self, target_gene_id):
       
        regulator_data = np.array([self.rna_data.get_values(REGULATOR_GENE, node_id).values for node_id in self.fate_map.node_leaves])
        target_data = np.array([self.rna_data.get_values(TARGET_GENE, node_id=node_id)[[target_gene_id]] for node_id in self.fate_map.node_leaves])
        return regulator_data, target_data

    def _reshape_atac_prior_input(self, target_gene_id):
        
        samples = [self.atac_data.get_node_data(node_id=node_id) for node_id in self.fate_map.node_leaves]
        mdata = np.abs([self.atac_data.reshape(sample, self.rna_data.regulator_genes, self.rna_data.target_genes)[[target_gene_id]].values for sample in samples])
        component1 = -np.log(1 + np.exp(-BETA_1 - BETA_2 * mdata))
        component2 = (-BETA_1 - BETA_2 * mdata) - np.log(1 + np.exp(-BETA_1 - BETA_2 * mdata))
        mdata = np.concatenate((component1, component2), axis=-1)
        return mdata

    def _reshape_lineage_prior_input(self):
      
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

    def infer_leaf_grns(self, target_gene_id: 'str' = None):
        """
        Infer GRN values for all leaf nodes for a specific target gene.

        **Parameters**

        :param str target_gene_id: Identifier of the target gene.

        **Returns**

        :return: Optimized GRN values after inference.
        :rtype: np.ndarray
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

    def infer_ancestor_grns(self, grn_values: 'Dict[str, List[float]]', target_gene_id: 'str'):
        """
        Infer ancestor GRN values by aggregating information from leaf nodes upward.

        **Parameters**

        :param Dict[str, List[float]] grn_values: Dictionary of GRN values for leaf nodes.
        :param str target_gene_id: Identifier of the target gene.

        **Returns**

        :return: Updated dictionary including ancestor GRN values.
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
        if not os.path.exists(self.saved_dir):
            os.makedirs(self.saved_dir)

        saved_path = os.path.join(self.saved_dir, f"target_gene_{target_gene_id}.csv")
        with open(saved_path, 'w+') as f:
            for node_id, grn_value in grn_values.items():
                line = [target_gene_id, node_id] + grn_value
                line = json.dumps(line)
                f.write(line + '\n')

    def infer_grns(self, max_processing):
        """
        Estimate GRN values for all target genes using parallel processing if specified.

        **Parameters**

        :param int max_processing: Maximum number of parallel processes to use.
        """
        if max_processing > 1:
            with Pool(max_processing) as p:
                p.map(self._get_strength, self.rna_data.target_genes)
        else:
            for target_gene_id in self.rna_data.target_genes:
                self._get_strength(target_gene_id)

    def get_target_networks(self, threshold: 'float', reverse=True):
        """
        Retrieve gene regulatory networks for each target gene based on a threshold.

        **Parameters**

        :param float threshold: Threshold for filtering GRN values.
        :param bool reverse: Flag to reverse direction if needed (unused in logic).

        **Returns**

        :return: Nested dictionary of GRN networks keyed by node and target gene.
        :rtype: Dict[str, Dict[str, Dict]]
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