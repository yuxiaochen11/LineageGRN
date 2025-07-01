API
=================


Inference
------------
.. list-table::
  :header-rows: 1

  * - **Main Function**
    - **Description**
  * - :doc:`construct_fate_map <construct_fate_map>`
    - Function for reconstructing time-scaled cell fate map.
  * - :doc:`GRNInference.infer_leaf_grn <GRNInference_module>`
    - A method of the class ```GRNInference``` to infer gene regulatory networks for the terminal cell types of the cell fate map.
  * - :doc:`GRNInference.infer_ancestor_grn <GRNInference_module>`
    - A method of the class ```GRNInference``` to infer gene regulatory networks for the progenitor cell types of the cell fate map.
  * - :doc:`GRNInference.infer_grn <GRNInference_module>`
    - A method of the class ```GRNInference``` to infer the dynamic gene regulatory networks for cell types of the cell fate map.
  * - :doc:`GRNInference.get_target_networks <GRNInference_module>`
    - A method of the class ```GRNInference``` to infer the target networks for cell types of a spesific target gene.


Analysis
------------
.. list-table::
   :header-rows: 1

   * - **Main Function**
     - **Description**
   * - :doc:`get_regulators_for_target_gene <get_regulators_for_target_gene>`
     - Function for retrieving the names and regulatory strengths of genes regulating a specified target gene.
   * - :doc:`get_targets_for_regulator_gene <get_targets_for_regulator_gene>`
     - Function for obtaining the target gene names and regulatory strengths regulated by a specified regulator gene.
   * - :doc:`identify_key_genes_differentiation <identify_key_genes_differentiation>`
     - Function for identifying key genes that drive the differentiation of specific progenitor cell types into downstream states.
   * - :doc:`identify_key_genes_fate_bias <identify_key_genes_fate_bias>`
     - Function for identifying key genes that drive the fate bias of progeny cells derived from a specific ancestral cell type.
   * - :doc:`cluster_regulatory_interactions <cluster_regulatory_interactions>`
     - Function for clustering regulatory interactions that coexist in similar cell types.
   * - :doc:`identify_regulatory_interactions_specificity <identify_regulatory_interactions_specificity>`
     - Function for identifying the specificity and constitutiveness of various clusters of regulatory interactions.




Plotting
----------
.. list-table::
   :header-rows: 1

   * - **Main Function**
     - **Description**
   * - :doc:`plot_regulatory_interactions_along_fatemap <plot_regulatory_interactions_along_fatemap>`
     - Function for plotting the dynamic changes in the number of regulatory interactions in the gene regulatory networks along the fate map.
   * - :doc:`plot_regulatory_genes_along_fatemap <plot_regulatory_genes_along_fatemap>`
     - Function for plotting the dynamic changes in the number of target genes regulated by a specific regulatory gene in the gene regulatory networks along the fate map.
   * - :doc:`plot_target_genes_along_fatemap <plot_target_genes_along_fatemap>`
     - Function for plotting the dynamic changes in the number of target genes regulated by a specific regulatory gene in the gene regulatory networks along the fate map.
   * - :doc:`plot_regulatory_network_along_fatemap <plot_regulatory_network_along_fatemap>`
     - Function for plotting the dynamic changes in the regulatory network of a specific regulatory gene along the fate map.
   * - :doc:`plot_regulatory_strength_along_fatemap <plot_regulatory_strength_to_target_gene_along_fatemap>`
     - Function for plotting the dynamic change in regulatory strength of each regulatory gene along a specific path of the fate map for a specific target gene.
   * - :doc:`plot_regulator_activity_across_lineages <plot_regulatory_gene_activity_across_lineages>`
     - Function for plotting a bar chart of the activity (number of regulated target genes) of a specific regulatory gene across different lineages and cell clusters.
   * - :doc:`plot_key_genes_differentiation <plot_key_genes_differentiation>`
     - Function for plotting the results of identifying key genes that drive the differentiation of specific progenitor cell types into downstream states.
   * - :doc:`plot_key_genes_fate_bias <plot_key_genes_fate_bias>`
     - Function for plotting the results of identifying key genes that drive the fate bias of progeny cells derived from a specific ancestral cell type.
   * - :doc:`plot_regulatory_interactions_clustering <plot_regulatory_interactions_clustering>`
     - Function for plotting the membership functions of the fuzzy clustering results of regulatory interactions, visualizing the degree of association between each regulatory interaction and its respective cluster.
   * - :doc:`plot_regulatory_interactions_in_celltypes <plot_regulatory_interactions_in_celltypes>`
     - Function for plotting the distribution of regulatory interactions within specific or constitutive regulatory clusters across different lineages or cell types.

