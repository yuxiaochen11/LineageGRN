def plot_regulatory_interactions_along_fatemap(grns_dict, path, output_path, figsize=(1.5, 1)):
    """
    Plot the number of dynamic edges for nodes along a specified path in the gene regulatory network.

    :param grns_dict: **A dictionary containing gene regulatory networks** for each node.
    :type grns_dict: dict
    :param path: **A list of node IDs** representing the path to trace.
    :type path: list of str
    :param output_path: **The path to save the output plot**.
    :type output_path: str
    :param figsize: **The size of the plot** (default is (1.5, 1)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the plot to the specified output path and displays it.
    :rtype: None

    """
    node_id_list = []   
    edge_number_list = [] 
    
    # Iterate through the nodes in the specified path
    for node_id in path:
        # Count the number of dynamic edges (negative regulatory interactions)
        edge_num = (grns_dict[node_id] < 0).sum().sum()
        node_id_list.append(node_id)       
        edge_number_list.append(edge_num)   
    
    # Create a DataFrame for the plot data
    edges_number_df = pd.DataFrame({
        'node_id': node_id_list, 
        'edge_number': edge_number_list
    })

    # Plot the number of dynamic edges
    plt.figure(figsize=figsize)  
    plt.plot(
        'node_id',
        'edge_number',
        data=edges_number_df,
        linestyle='-',
        marker='o',
        linewidth=0.5,
        markersize=3,
        color='#2EA7E0',          
        markerfacecolor='none'    
    )
    
    # Label the axes with appropriate font size
    plt.xlabel('Cell type', fontsize=6, fontname='Arial')   
    plt.ylabel('# Regulatory interaction', fontsize=6, fontname='Arial')  
    
    # Set tick sizes and rotation for better readability
    plt.xticks(rotation=-90, fontsize=5, fontname='Arial')
    plt.yticks(fontsize=5, fontname='Arial')
    
    # Save the plot in EPS format and display it
    plt.savefig(output_path + 'dynamic_edge_number.eps', format='eps', bbox_inches='tight')
    plt.show()

def plot_regulatory_genes_along_fatemap(input_path, regulator_names, regulator_id, path, threshold, output_path, figsize=(1.5, 1)):
    """
    Plot the number of target genes regulated by a specific regulator along a specified path.

    :param input_path: **The path to the input data files**.
    :type input_path: str
    :param regulator_names: **A list of regulator names**.
    :type regulator_names: list of str
    :param regulator_id: **The ID of the regulator to analyze**.
    :type regulator_id: str
    :param path: **A list of node IDs representing the path to trace**.
    :type path: list of str
    :param threshold: **A threshold for determining the significance of target gene interactions**.
    :type threshold: float
    :param output_path: **The directory where the plot will be saved**.
    :type output_path: str
    :param figsize: **The size of the plot** (default is (1.5, 1)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the plot to the specified output path and displays it.
    :rtype: None

    """
    
    node_id_list = []      
    target_number_list = []  

    # Iterate through the nodes in the specified path
    for node_id in path:
        # Retrieve the target genes regulated by the specified regulator
        target_dict = get_target_genes(regulator_id, node_id, input_path, regulator_names)
        
        # Count the number of target genes that meet the threshold
        target_num = sum(1 for v in target_dict.values() if abs(v) > threshold)
        
        # Append the results to the lists
        node_id_list.append(node_id)
        target_number_list.append(target_num)

    # Create a DataFrame for the plot data
    target_number_df = pd.DataFrame({'node_id': node_id_list, 'target_number': target_number_list})

    # Create and plot the figure
    plt.figure(figsize=figsize) 
    plt.plot(
        'node_id', 
        'target_number', 
        data=target_number_df, 
        linestyle='-', 
        marker='o', 
        linewidth=0.5, 
        markersize=3, 
        color='#2EA7E0',           
        markerfacecolor='none'     
    )
    
    # Label the axes with appropriate font size
    plt.xlabel('Cell type', fontsize=6, fontname='Arial')
    plt.ylabel('# Target gene', fontsize=6, fontname='Arial')
    
    # Set tick sizes and rotation for better readability
    plt.xticks(rotation=-90, fontsize=5, fontname='Arial')
    plt.yticks(fontsize=5, fontname='Arial')

    # Save the plot as EPS and display it
    plt.savefig(output_path + f'dynamic_target_gene_number_for_{regulator_id}.eps', format='eps', bbox_inches='tight')
    plt.show()

def plot_target_genes_along_fatemap(target_gene_names, regulator_names, regulatory_module, input_path, nodes, threshold, output_path, figsize=(2.2, 1.5)):
    """
    Plot a heatmap showing the number of regulators for each target gene across specified nodes.

    :param target_gene_names: **A list of target gene IDs to analyze**.
    :type target_gene_names: list of str
    :param regulator_names: **A list of regulator names to be used for analysis**.
    :type regulator_names: list of str
    :param regulatory_module: **The regulatory module to analyze**. Options are `'positive'`, `'negative'`, or `'total'`.
    :type regulatory_module: str
    :param input_path: **The folder path containing input files with regulatory data**.
    :type input_path: str
    :param nodes: **A list of node IDs for which to calculate the number of regulators**.
    :type nodes: list of str
    :param threshold: **The threshold value for considering a regulator active**.
    :type threshold: float
    :param output_path: **The path where the output plot will be saved**.
    :type output_path: str
    :param figsize: **The size of the plot** (default is (2.2, 1.5)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the heatmap to the specified output path and displays it.
    :rtype: None

    """
    
    plt.rcParams['font.size'] = 6  # Set font size for the plot

    regulator_number_df = pd.DataFrame()  # DataFrame to hold the regulator counts for each target gene
    
    # Iterate through each target gene
    for target_gene_id in target_gene_names:
        regulator_number_list = []  # List to store regulator counts for the current target gene
        
        # Iterate through each node
        for node_id in nodes:
            # Get the regulators for the target gene at the current node
            regulator_dict = get_regulators(target_gene_id, node_id, input_path, regulator_names)
            
            # Count the number of regulators based on the regulatory module (positive, negative, or total)
            if regulatory_module == 'positive':
                regulator_num = sum(1 for v in regulator_dict.values() if v > threshold)
                color = 'Reds'  # Set color map for positive regulation
            elif regulatory_module == 'negative':
                regulator_num = sum(1 for v in regulator_dict.values() if v < -threshold)
                color = 'Blues'  # Set color map for negative regulation
            else:
                regulator_num = sum(1 for v in regulator_dict.values() if abs(v) > threshold)
                color = 'Oranges'  # Set color map for total regulation
            
            regulator_number_list.append(regulator_num)  # Append the count to the list
        
        # Add the regulator count list for the current target gene to the DataFrame
        regulator_number_df_for_target_id = pd.DataFrame({target_gene_id: regulator_number_list})
        regulator_number_df = pd.concat([regulator_number_df, regulator_number_df_for_target_id], axis=1)
    
    regulator_number_df.index = nodes  # Set the index of the DataFrame to the node IDs
    
    # Create the heatmap plot
    plt.figure(figsize=figsize)
    g = sns.heatmap(data=regulator_number_df.T, cmap=plt.get_cmap(color))  # Transpose DataFrame for heatmap
    
    # Adjust tick labels and font size
    plt.setp(g.get_yticklabels(), fontsize=5, fontname='Arial')
    plt.setp(g.get_xticklabels(), fontsize=5, fontname='Arial', rotation=-90)
    g.tick_params(axis='both', which='major', labelsize=5, length=1)
    
    # Set colorbar ticks size
    if g.collections:
        g.collections[0].colorbar.ax.tick_params(labelsize=5, length=1)
    
    # Label the axes
    plt.xlabel('Cell type', fontsize=6, fontname='Arial')
    plt.ylabel('Target genes', fontsize=6, fontname='Arial')

    # Save the plot as EPS and display it
    plt.savefig(output_path + 'dynamic_regulator_number_heatmap.eps', format='eps', bbox_inches='tight')
    plt.show()

def plot_regulatory_strength_to_target_gene_along_fatemap(input_path, path, regulator_names, target_gene_id, output_path, figsize=(3, 1.5)):
    """
    Plot the regulatory strength of different regulators for a specified target gene across various nodes.

    :param input_path: **The folder path containing input files with regulatory data**.
    :type input_path: str
    :param path: **A list of node IDs to analyze**.
    :type path: list of str
    :param regulator_names: **A list of regulator names to be used as column headers**.
    :type regulator_names: list of str
    :param target_gene_id: **The ID of the target gene for which to plot regulatory strengths**.
    :type target_gene_id: str
    :param output_path: **The path where the output plot will be saved**.
    :type output_path: str
    :param figsize: **The size of the plot** (default is (3, 1.5)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the plot to the specified output path and displays it.
    :rtype: None

    """

    file_path =  input_path + '/target_gene_' + target_gene_id + '.csv'
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
            
        data = [ast.literal_eval(line.strip()) for line in lines]
        target_gene_data = pd.DataFrame(data).iloc[:, 3:]
        target_gene_data.columns = regulator_names
        target_gene_data.index = list(pd.DataFrame(data).iloc[:, 1])

    scatter_df = []  

    # Loop over each regulator and node to extract the regulatory strength data
    for i in range(target_gene_data.shape[1]):
        for j in range(target_gene_data.shape[0]):
            row = [target_gene_data.columns[i], abs(target_gene_data.iloc[j, i]), target_gene_data.index[j]]
            scatter_df.append(row)
    
    scatter_df = pd.DataFrame(scatter_df, columns=['regulator', 'value', 'node_id'])
    scatter_df = scatter_df[scatter_df['node_id'].isin(path)]  # Filter for specified nodes
    
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("#FFFFFF")  
    ax.set_facecolor("#FFFFFF") 
    ax.axhline(0.3, color='#231815', ls=(0, (5, 5)), alpha=1, zorder=0,linewidth=0.5)  
    nodes = sorted(scatter_df["node_id"].unique())

    # Plot the regulatory strength for each node with a scatter plot
    for node, color, edgecolors, marker in zip(nodes, 
                                    ['#FAD5D2','#C7E8FA', '#FCFAD3',  '#D2E8CA', '#2ca02c', '#D0F0C0', '#FCFAD3'],
                                    ['#E60012','#2171A9', '#EF7C20',  '#339939', '#1f77b4', '#98FB98', '#FCFAD3'],
                                    ["D",'o',  "s", "*",'h','x']):
        data = scatter_df[scatter_df["node_id"] == node]
        ax.scatter("regulator", "value", s=14, color=color, edgecolors=edgecolors, linewidths=0.6, marker=marker, alpha=1, data=data)

    ax.tick_params(axis='x', rotation=-90)  # Rotate x-axis tick labels

    fig.suptitle(target_gene_id, x=0.5, y=0.975, ha="center", fontsize=8, color='#231815')

    # Add legend for nodes
    legend = ax.legend(loc=(0.05, 1.05), labelspacing=0.2, markerscale=0.8, frameon=False)
    for text, node in zip(legend.get_texts(), nodes):
        text.set_text(node)
        text.set_fontsize(6)
        text.set_fontname('Arial')

    legend.set_title("Cell type")
    legend_title = legend.get_title()
    legend_title.set_fontsize(6.5)
    legend_title.set_fontname('Arial')
    legend_title.set_ha("left")
    
    # Set plot appearance and formatting
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["left"].set_color('#231815')
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_color('#231815')
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(length=0.2) 
    
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0], size=5, fontname='Arial')
    ax.set_ylabel("Regulatory strength", size=6, fontname='Arial')
    
    xticks = list(scatter_df['regulator'].unique())
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, size=5, fontname='Arial')
    ax.set_xlabel("Regulatory gene", size=6, fontname='Arial')

    # Save and show the plot
    plt.savefig(output_path + 'dynamic_regulatory_strength.eps', format='eps', bbox_inches='tight')
    plt.show()

def plot_regulatory_network_along_fatemap(regulator_id, grns_dict, nodes, output_path, threshold=0.1, figsize=(6.8, 1.8)):
    """
    Plot the dynamic regulatory network for a specific regulator across multiple nodes.

    :param regulator_id: **The ID of the regulator to analyze**.
    :type regulator_id: str
    :param grns_dict: **A dictionary containing gene regulatory networks for each node**.
    :type grns_dict: dict
    :param nodes: **A list of node IDs for which the regulatory network is to be plotted**.
    :type nodes: list of str
    :param output_path: **The path to save the output plot**.
    :type output_path: str
    :param threshold: **The threshold for considering a regulatory interaction as significant** (default is 0.1).
    :type threshold: float, optional
    :param figsize: **The size of the plot** (default is (6.8, 1.8)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the plot to the specified output path and displays it.
    :rtype: None
    """

    def custom_layout(G, edge_lengths, radius=1.0):
        pos = {}
        nodes_list = list(G.nodes())
        num_regulators = len(nodes_list) - 1
        if num_regulators == 0:
            return pos
        angle_step = 2 * np.pi / num_regulators
        pos[nodes_list[0]] = np.array([0, 0])
        for i, rid in enumerate(nodes_list[1:]):
            angle = i * angle_step
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pos[rid] = np.array([x, y])

        for (u, v), length in zip(G.edges(), edge_lengths):
            if u == nodes_list[0] or v == nodes_list[0]:
                if u == nodes_list[0]:
                    rid = v
                else:
                    rid = u
                norm_length = length / max(edge_lengths)
                pos[rid] = np.array([
                    radius * norm_length * np.cos(angle_step * (nodes_list[1:].index(rid))),
                    radius * norm_length * np.sin(angle_step * (nodes_list[1:].index(rid)))
                ])
        return pos

    def color_map(value):
        if value > 0:
            return plt.cm.Reds(value * 1)
        else:
            return plt.cm.Blues(-value * 2)

    def merge_and_extract_regulator(nodes_list, grns_dictionary, rid):
        df_list = [grns_dictionary[node] for node in nodes_list]
        merged_df = pd.concat(df_list, axis=0)
        return merged_df.loc[:, rid]

    def edge_colormap_plot(ax, regulator_id, grns_dictionary, node, threshold, global_min, global_max, ticks):
        value_df = pd.DataFrame(grns_dictionary[node].loc[:, regulator_id])
        value_all = merge_and_extract_regulator(nodes, grns_dictionary, regulator_id)

        targets_all = list(set(value_all[abs(value_all) > threshold].index))
        targets = list(value_df[abs(value_df[regulator_id]) > threshold].index)  # 修正此处
        values = list(value_df.loc[abs(value_df[regulator_id]) > threshold, regulator_id])

        data_all = pd.DataFrame({'regulator': [regulator_id] * len(targets_all), 'targets': targets_all})
        data = {'regulator': [regulator_id] * len(targets), 'targets': targets, 'value': values}
        data_df = pd.DataFrame(data)

        graph_df = pd.merge(data_all, data_df, on=['regulator', 'targets'], how='left')
        graph_df['value'] = graph_df['value'].fillna(0)
        graph_df = graph_df.loc[graph_df['targets'] != regulator_id]

        G = nx.from_pandas_edgelist(graph_df, 'regulator', 'targets')
        edge_lengths = [0.1] * len(graph_df['value'])
        pos = custom_layout(G, edge_lengths, radius=10)

        edge_colors = [color_map(c * 1) for c in graph_df['value']]
        node_colors = ['#E4F0CF' if target == regulator_id else '#DCDCDD' for target in G.nodes()]
        node_sizes = [200 if target == regulator_id else 100 for target in G.nodes()]

        # Draw the network
        nx.draw(
            G, pos, node_color=node_colors, edge_color=edge_colors,
            width=1, with_labels=True, node_size=node_sizes,
            font_size=4, ax=ax
        )

        norm = Normalize(vmin=global_min, vmax=global_max)
        color_list = [color_map(v) for v in np.linspace(global_min, global_max, 200)]
        cmap = ListedColormap(color_list)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Regulatory strength', fontsize=6, fontname="Arial")
        cbar.set_ticks(ticks)
        cbar.set_ticklabels([f'{tick:.2f}' for tick in ticks], fontname="Arial")
        cbar.ax.tick_params(labelsize=5)

        ax.set_title(f'{node}', fontsize=7, fontname="Arial")

    # Create subplots for each node
    fig, axs = plt.subplots(1, len(nodes), figsize=figsize)
    axs = axs.flatten()
    global_min = float('inf')
    global_max = float('-inf')

    # Find global min and max for color normalization
    for node in nodes:
        values = list(grns_dict[node].loc[:, regulator_id])
        current_min = min(values)
        current_max = max(values)
        global_min = min(global_min, current_min)
        global_max = max(global_max, current_max)

    ticks = np.linspace(global_min, global_max, num=6)

    # Plot for each node
    for i, node in enumerate(nodes):
        if i < len(axs):
            edge_colormap_plot(axs[i], regulator_id, grns_dict, node,
                               threshold, global_min, global_max, ticks)

    plt.tight_layout()
    fig.savefig(output_path + 'dynamic_regulatory_network.eps', format='eps', bbox_inches='tight')
    plt.show()

def plot_regulatory_gene_activity_across_lineages(edges_dict, input_folder_path, regulator_id, regulator_names, threshold, output_path, figsize=(3, 3)):
    """
    Plot the activity of a specified regulator across different lineages in a polar plot.

    :param edges_dict: **A dictionary containing edges for different lineages**.
    :type edges_dict: dict
    :param input_folder_path: **The folder path where input data is located**.
    :type input_folder_path: str
    :param regulator_id: **The ID of the regulator to analyze**.
    :type regulator_id: str
    :param regulator_names: **A list of regulator names to analyze**.
    :type regulator_names: list of str
    :param threshold: **The threshold for considering a regulatory interaction as significant**.
    :type threshold: float
    :param output_path: **The path to save the output plot**.
    :type output_path: str
    :param figsize: **The size of the polar plot** (default is (3, 3)).
    :type figsize: tuple of float, optional

    :returns: **None**: This function saves the plot to the specified output path and displays it.
    :rtype: None

    """

    fate_map_dict = {}
    for lineage in edges_dict.keys():
        edges = parse_edge_dict(edges_dict[lineage])  # Parse edges for each lineage
        fate_map = FateMap(edges)  # Create a FateMap object for the lineage
        fate_map_dict.update({lineage: fate_map})

    lineage_df = pd.DataFrame()

    # Process each lineage to calculate the number of targets for each regulator
    for lineage in edges_dict.keys():
        input_path = os.path.join(input_folder_path, lineage)
        nodes = fate_map_dict[lineage].nodes.keys()
        for node_id in nodes:
            target_number_list = [lineage, node_id]
            for reg_id in regulator_names:
                target_dict = get_target_genes(reg_id, node_id, input_path, regulator_names)
                target_num = sum(1 for v in target_dict.values() if abs(v) > threshold)
                target_number_list.append(target_num)
            data = pd.DataFrame({node_id: target_number_list})
            lineage_df = pd.concat([lineage_df, data], axis=1)

    lineage_df.index = ['lineage', 'node_id'] + regulator_names
    lineage_df = lineage_df.T

    def add_labels(angles, values, labels, offset, ax):
        """
        Add labels to the polar plot at specified angles and values.

        :param angles: **Angles at which to place the labels**.
        :type angles: np.ndarray
        :param values: **Values at which to place the labels**.
        :type values: np.ndarray
        :param labels: **Labels to display**.
        :type labels: np.ndarray
        :param offset: **Angle offset for positioning labels**.
        :type offset: float
        :param ax: **The axes to add labels to**.
        :type ax: matplotlib.axes.Axes

        :returns: **None**: This function adds labels directly to the axes.
        :rtype: None
        """
        padding = 1  # Padding for label positioning
        for angle, value, label in zip(angles, values, labels):
            rotation = np.rad2deg(angle + offset)
            if angle <= np.pi:
                alignment = "right"
                rotation += 180
            else:
                alignment = "left"

            ax.text(
                x=angle, 
                y=value + padding, 
                s=label, 
                ha=alignment, 
                va="center", 
                rotation=rotation, 
                rotation_mode="anchor",
                fontsize=5,
                fontname='Arial'
            )

    # Extract the values, labels, and groups for plotting
    VALUES = lineage_df[regulator_id].values
    LABELS = lineage_df["node_id"].values
    GROUP = lineage_df["lineage"].values
    OFFSET = np.pi / 2  # Angle offset for the plot
    PAD = 2  # Padding between groups
    ANGLES_N = len(VALUES) + PAD * len(np.unique(GROUP))
    ANGLES = np.linspace(0, 2 * np.pi, num=ANGLES_N, endpoint=False)
    WIDTH = (2 * np.pi) / len(ANGLES)  # Width of the bars
    GROUPS_SIZE = [len(i[1]) for i in lineage_df.groupby("lineage")]

    offset = 0
    IDXS = []

    # Indexing for plotting bars in the correct order
    for size in GROUPS_SIZE:
        IDXS += list(range(offset + PAD, offset + size + PAD))
        offset += size + PAD

    fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})
    ax.set_theta_offset(OFFSET)  # Set angle offset
    ax.set_ylim(-20, 40)  # Set the radial limits
    ax.set_frame_on(False)  # Remove the plot frame
    ax.xaxis.grid(False)  # Hide x-axis grid
    ax.yaxis.grid(False)  # Hide y-axis grid
    ax.set_xticks([])  # Hide x-ticks
    ax.set_yticks([])  # Hide y-ticks

    COLORS_GROUP = ['#48ABE1', '#B28247', '#CA87B8', '#EC6655', '#78BF5B', '#F6B956']
    COLORS = [COLORS_GROUP[i] for i, size in enumerate(GROUPS_SIZE) for _ in range(size)]

    # Plot bars representing the regulator activity
    ax.bar(
        ANGLES[IDXS], VALUES, width=WIDTH, color=COLORS, 
        edgecolor="white", linewidth=0.5
    )

    # Add labels to the plot
    add_labels(ANGLES[IDXS], VALUES, LABELS, OFFSET, ax)

    # Draw lineages separation and add lineage labels
    offset = 0
    for lineage, size in zip(edges_dict.keys(), GROUPS_SIZE):
        x1 = np.linspace(ANGLES[offset + PAD], ANGLES[offset + size + PAD - 1], num=50)
        ax.plot(x1, [-5] * 50, color="#E1E5E5")

        ax.text(
            np.mean(x1), -10, lineage, color="#333333", fontsize=5, ha="center", va="center", fontname="Arial"
        )

        x2 = np.linspace(ANGLES[offset], ANGLES[offset + PAD - 1], num=50)
        for y_val in [5, 10, 15, 20, 25, 30]:
            ax.plot(x2, [y_val] * 50, color="#bebebe", lw=0.1)
        
        offset += size + PAD

    fig.show()  # Display the plot
    fig.savefig(output_path + 'dynamic_regulator_activity.eps', format='eps', bbox_inches='tight')


def plot_key_genes_differentiation(data, nodes, regulator_names, output_path, figsize=(5, 1)):
    """
    Plots key regulators based on their regulatory strength and node information.

    :param data: **DataFrame containing regulatory information including probabilities and activation scales**.
    :type data: pandas.DataFrame
    :param nodes: **List of node IDs to filter the data**.
    :type nodes: list of str
    :param regulator_names: **List of regulator names for labeling the plot**.
    :type regulator_names: list of str
    :param output_path: **Path where the output plot will be saved**.
    :type output_path: str
    :param figsize: **The size of the plot** (default is (5, 1)).
    :type figsize: tuple of float, optional

    :returns: **ax**: The matplotlib axis object with the plot.
    :rtype: matplotlib.axes.Axes

    """

    data = data[data['node_id'].isin(nodes)]

    fig, ax = plt.subplots(figsize=figsize)

    # Define colors and normalization for the plot
    COLORS = ['#FFFDDF', "#7FCDBB", "#225EA8", '#132860']
    cmap = mcolors.LinearSegmentedColormap.from_list("colormap", COLORS, N=100)
    norm = plt.Normalize(vmin=0.2, vmax=data['DRS'].max())

    # Define bubble sizes for legend
    size_values = [1, 0.75, 0.5]  
    size_labels = ['1', '0.75', '0.5'] 
    size_colors = ['black', 'black', 'black'] 
    size_edgecolors = ['black', 'black', 'black'] 

    # Plot the legend for bubble sizes
    for i, size in enumerate(size_values):
        if size <= 0.5:
            size = size * 0.25
        if (size <= 0.75 and size > 0.5):
            size = size * 0.5
        ax.scatter([i], [-2], s=size * 30, color=size_colors[i], edgecolor=size_edgecolors[i], label=size_labels[i])

    # Plot regulators as scatter points
    for i, regulator in enumerate(regulator_names):
        d = data[data['regulator_id'] == regulator]
        y = d['node_id']
        x = [i] * len(y)  # Each regulator has its own x-position
        color = cmap(norm(d["DRS"]))  # Color scale based on DRS
        sizes_data = d['DRC'] * 30  # Adjust size based on DRC
        sizes_data = [x * 0.5 if x <= 0.75 else x for x in sizes_data]  
        sizes_data = [x * 0.5 if x < 0.5 else x for x in sizes_data] 

        ax.scatter(x, y, color=color, s=sizes_data, edgecolor='black', linewidths=0.25)

    # Configure plot appearance
    ax.set_frame_on(False)
    ax.grid(linewidth=0.5, alpha=0.7, color='#E5E5E5')
    ax.set_axisbelow(True)
    ax.set_xticks(np.arange(len(regulator_names)))
    ax.set_xticklabels(regulator_names, rotation=90, color='black', fontsize=5, fontname='Arial', ha='center')
    ax.tick_params(axis='x', which='both', length=0)  

    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, color='black', fontsize=5, fontname='Arial', va='center')
    ax.tick_params(axis='y', which='both', length=0)  

    ax.set_xlabel('Regulatory gene', color='black', fontsize=6, fontname='Arial')
    ax.set_ylabel("Cell type", color='black', fontsize=6, fontname='Arial')

    # Shrink y-limits to avoid labels being cut off
    y_shrunk = 0.01
    y_lower, y_upper = ax.get_ylim()
    ax.set_ylim(y_lower + y_shrunk, y_upper - y_shrunk)

    # Color bar for DRS values
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar_ax = fig.add_axes([1, 0.45, 0.02, 0.4])  # Position of the colorbar
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label('DRS', fontsize=6, fontname='Arial')
    cbar.ax.tick_params(labelsize=5, length=0.5)  

    # Add legend for DRC sizes
    ax.legend(loc='upper right', bbox_to_anchor=(1.7, 1), fontsize=5, frameon=False, title='DRC', title_fontsize=6)

    # Save and show the plot
    fig.savefig(output_path + 'find_diff_genes.eps', format='eps', bbox_inches='tight')
    plt.show()
    
    return ax


def plot_key_genes_fate_bias(df, child_nodes, output_path, figsize=(1.5, 12)):
    """
    Plots the regulatory gene bias for two child nodes based on their regulatory strength and bias values.

    This function creates a scatter plot showing the bias between two child nodes (`child_nodes[0]` and `child_nodes[1]`),
    using regulatory gene bias (`FBP`) and regulatory strength (`FBI`) values to position and color the data points.
    The size of the points reflects the `FBP` values, and a colorbar is added to represent the regulatory strength (`BRS`).

    :param df: **DataFrame containing the regulatory data.**
    :type df: pandas.DataFrame
    :param child_nodes: **List of two child nodes to compare.**
    :type child_nodes: list of str
    :param output_path: **Path where the output plot will be saved.**
    :type output_path: str
    :param figsize: **Size of the figure.** Default is (1.5, 12).
    :type figsize: tuple of float, optional

    :returns: None
    :rtype: None

    """

    # Extract relevant columns for the two child nodes
    FBP_0 = 'FBP_' + child_nodes[0]
    FBP_1 = 'FBP_' + child_nodes[1]
    FBI_0 = 'FBI_' + child_nodes[0]
    FBI_1 = 'FBI_' + child_nodes[1]

    # Filter out rows where both biases are less than 1
    df = df[~((df[FBP_0] < 1) & (df[FBP_1] < 1))].copy()
    df.reset_index(drop=True, inplace=True)

    # Prepare data for plotting
    genes = df['Regulatory gene'].values
    y_vals = np.arange(len(genes))  # Y positions for each gene
    x_endo, x_epi = 0.0, 0.3  # X positions for two child nodes
    positions = [0.0, 0.4999, 0.5000, 0.5001, 1.0]  # Color transition positions
    colors = [
        "#FFFFFF",  # pos=0.0
        "#7FCDBB",  # pos=0.4999
        "#7FCDBB",  # pos=0.5000
        "#225EA8",  # pos=0.5001
        "#132860",  # pos=1.0
    ]
    my_cmap = LinearSegmentedColormap.from_list("my_cmap", list(zip(positions, colors)))

    # Determine global minimum and maximum values for both FBP and FBI
    FBI_min = df[[FBI_0, FBI_1]].min().min()
    FBI_max = df[[FBI_0, FBI_1]].max().max()
    FBP_global_min = df[[FBP_0, FBP_1]].min().min()
    FBP_global_max = df[[FBP_0, FBP_1]].max().max()

    # Function to convert FBP values into marker sizes
    def FBP_to_size(FBP_value):
        size_min_small = 1
        size_max_small = 10
        size_max_large = 40
        if FBP_value < 1:
            return size_min_small + (size_max_small - size_min_small) * (FBP_value - 0) / (1 - 0)
        else:
            if FBP_global_max == 1:
                return size_max_small
            else:
                return size_max_small + (size_max_large - size_max_small) * \
                       (FBP_value - 1) / (FBP_global_max - 1)
    
    # Apply size function to the data
    sizes_endo = df[FBP_0].apply(FBP_to_size)
    sizes_epi = df[FBP_1].apply(FBP_to_size)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize)  
    sc_endo = ax.scatter(
        np.full_like(y_vals, x_endo, dtype=float),  
        y_vals,
        c=df[FBI_0],
        s=sizes_endo,
        cmap=my_cmap,
        vmin=FBI_min,
        vmax=FBI_max,
        edgecolors='k',
        alpha=1,
        linewidths=0.25
    )
    sc_epi = ax.scatter(
        np.full_like(y_vals, x_epi, dtype=float),   
        y_vals,
        c=df[FBI_1],
        s=sizes_epi,
        cmap=my_cmap,
        vmin=FBI_min,
        vmax=FBI_max,
        edgecolors='k',
        alpha=1,
        linewidths=0.25
    )

    # Customize axes and labels
    ax.set_xticks([x_endo, x_epi])
    ax.set_xticklabels(child_nodes, fontsize=5, fontname='Arial', rotation=0)
    ax.set_yticks(y_vals)
    ax.set_yticklabels(genes, fontsize=5, fontname='Arial')
    ax.set_ylim(-0.7, len(genes) - 0.7)
    ax.set_xlim(-0.2, x_epi + 0.2)

    # Customize plot appearance
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    ax.tick_params(width=0.25, length=1)

    # Add colorbar for BRS values
    cbar = plt.colorbar(sc_endo, ax=ax)
    cbar.ax.tick_params(labelsize=5, width=0.2)
    cbar.set_label("BRS", fontsize=5, fontname='Arial')
    cbar_ticks = [FBI_min, 0.5, FBI_max]
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{FBI_min:.2g}", "0.5", f"{FBI_max:.2g}"])

    # Add legend for FBP values
    FBP_legend_values = [0.5, 1.0, 2.0, FBP_global_max]  
    legend_handles = []
    legend_labels = []
    for val in FBP_legend_values:
        val_clamped = max(FBP_global_min, min(val, FBP_global_max))
        size_ = FBP_to_size(val_clamped)
        h = ax.scatter([], [], s=size_, c='gray', edgecolors='k', linewidths=0.25)
        legend_handles.append(h)
        legend_labels.append(f"{val_clamped:.2g}")
    
    FBP_legend = ax.legend(
        handles=legend_handles,
        labels=legend_labels[0:(len(legend_labels)-1)],
        title="BRN",
        bbox_to_anchor=(1.2, 1.0),
        title_fontsize=6,  
        prop={'family': 'Arial', 'size': 5}
    )
    ax.add_artist(FBP_legend)

    # Set labels
    ax.set_ylabel("Regulatory gene", fontsize=6, fontname='Arial')

    # Finalize the plot layout and save the output
    plt.tight_layout()
    plt.savefig(output_path + 'find_fate_bias_genes.eps', format='eps', bbox_inches='tight')
    plt.show()



def plot_regulatory_interactions_clustering(weight_matrix, output_path, width=3, height=7):
    """
    Plots a clustered heatmap of the weights in a given edge cluster.

    This function generates a heatmap representing the weights of edges in a given edge cluster.
    The heatmap is clustered using hierarchical clustering to group similar weights together.

    :param weight_matrix: **A matrix containing the weights of edges in the edge cluster.**
    :type weight_matrix: pandas.DataFrame
    :param output_path: **Path where the output plot will be saved.**
    :type output_path: str
    :param width: **Width of the output plot.** Default is 3.
    :type width: float, optional
    :param height: **Height of the output plot.** Default is 7.
    :type height: float, optional

    :returns: None
    :rtype: None

    """
    
    COLORS = ['#FFFDDF', "#7FCDBB", "#225EA8"]

    g = sns.clustermap(
        weight_matrix,
        cmap=mcolors.LinearSegmentedColormap.from_list("colormap", COLORS, N=100)
    )
    
    g.fig.set_size_inches(width, height)  
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=5, fontname='Arial')
    plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=5, fontname='Arial', rotation=50)

    g.ax_heatmap.tick_params(axis='both', which='major', labelsize=5, length=0)
    g.cax.tick_params(labelsize=5, length=0.5)

    plt.savefig(output_path + 'edge_cluster_weight', format='eps', bbox_inches='tight')

def plot_regulatory_interactions_in_celltypes(locate_edge_cluster_to_nodes_df, output_path, figsize=(1.8, 2.2)):
    """
    Plots a heatmap showing the association between edge clusters and nodes.

    This function generates a heatmap that visualizes the relationship between edge clusters and nodes.
    The heatmap's cells are annotated with the values from the DataFrame.

    :param locate_edge_cluster_to_nodes_df: **A DataFrame containing the relationship between edge clusters and nodes.**
    :type locate_edge_cluster_to_nodes_df: pandas.DataFrame
    :param output_path: **Path where the output plot will be saved.**
    :type output_path: str
    :param figsize: **Size of the figure.** Default is (1.8, 2.2).
    :type figsize: tuple of float, optional

    :returns: None
    :rtype: None

    """

    COLORS = ['#FFFDDF', "#7FCDBB", "#225EA8"]

    plt.figure(figsize=figsize)

    g = sns.heatmap(locate_edge_cluster_to_nodes_df, fmt=".3f", annot=True, annot_kws={"fontsize": 5, "fontname": "Arial"},
                    cmap=mcolors.LinearSegmentedColormap.from_list("colormap", COLORS, N=100))

    plt.setp(g.get_yticklabels(), fontsize=5, fontname='Arial')

    plt.setp(g.get_xticklabels(), fontsize=5, fontname='Arial', rotation=50)

    g.tick_params(axis='both', which='major', labelsize=5, length=0)

    if g.collections:
        g.collections[0].colorbar.ax.tick_params(labelsize=5, length=0.5)

    plt.xlabel('Edge clusters', fontsize=6, fontname='Arial')
    plt.ylabel('Cell type', fontsize=6, fontname='Arial')

    plt.savefig(output_path + 'edge_cluster_to_nodes', format='eps', bbox_inches='tight')