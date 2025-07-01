def construct_fate_map(sc_mat, cell_types, t_S, beta=1, max_iter=1000, tol=1e-5, output_edges=False):
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

