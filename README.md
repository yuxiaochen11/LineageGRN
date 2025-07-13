# LineageGRN
***
LineageGRN aims to infer dynamic gene regulatory networks at a real-time scale using lineage tracing, scRNA-seq, and scATAC-seq datasets. The inferred GRNs are dependency and linked by a cell fate map representing the developmental trajectory from a progenitor cell to their descendants (sampled cells), which is modeled by a Bayesian network. 
***

### Features
The currently available features include:
- Inferring time-scaled cell fate maps
- Inferring dynamic gene regulatory networks
- Revealing how inferred GRNs reconfigure along cell lineages
- Identifying key regulatory genes driving cell differentiation
- Identifying key regulatory genes steering cell fate bias
- Identifying specific and constitutive regulatory interactions.  


### Environment
python3.11

#### Installation
```shell
$ python setup.py install
```

***
### Usage overview

#### Inferring the real-time-scaled cell fate map
You have two options of inputs to infer the cell fate map. First, input lineage tracing data as a numpy array:

```python
from cell_lineage_reconstruction import construct_fate_map
import numpy as np

# Input the lineage tracing data
elements = np.array([0, 1, 2, 3, 4, 5])
sc_mat = np.random.choice(elements, size=(30, 15)) 
t_S=15

# Input the the cell type of each sample cell if you want each node of the cell fate map to be defined as a cell cluster with the same cell type.
cell_types = {
            "Cell0":"C","Cell1":"A","Cell2":"A","Cell3":"B","Cell4":"C",
            "Cell5":"A","Cell6":"B","Cell7":"D","Cell8":"D","Cell9":"B",
            "Cell10":"B","Cell11":"C","Cell12":"C","Cell13":"A","Cell14":"A",
            "Cell15":"A","Cell16":"D","Cell17":"D","Cell18":"D","Cell19":"B",
            "Cell20":"C","Cell21":"A","Cell22":"C","Cell23":"C","Cell24":"B",
            "Cell25":"B","Cell26":"A","Cell27":"A","Cell28":"D","Cell29":"B"
            }

#Reconstruct the real-time-scaled cell state tree as a fate map
fate_map = construct_fate_map(sc_mat, cell_types, t_S, beta=1, max_iter=1000, tol=1e-5)
```
You can also directly input the topology structure of the cell fate map on which you want to observe the dynamic changes of gene regulatory networks.

```python
from cell_lineage_reconstruction import FateMap

# Input the cell fate map directly
edge_dict = {"grn0->grn1":0.4,"grn0->grn2":0.1,"grn1->grn3":0.6,"grn1->grn4":0.6,"grn2->grn5":0.9,"grn2->grn6":0.9}
parse_edge = parse_edge_dict(edge_dict)
fate_map = FateMap(parse_edge)
```

#### Inferring the dynamic gene regulatory network
Input the real-time-scaled fate map inferred from the lineage tracing data and prepare your transcriptome and chromatin accessibility data in the format of the sample files, and you can start inferring the real-time resolved dynamic regulatory networks right away:

```python
from gene_regulatory_network import GRNInference
import pandas as pd

# Input the transcriptome and chromatin accessibility data
atac_file_path = "examples/data/simulation/input/atac_data.csv"
expression_file_path = "examples/data/simulation/input/expression_data.csv"

# Define an output path
saved_dir = 'examples/results/simulation/inferred_grn'

grn_inference = GRNInference(atac_file_path, expression_file_path,fate_map, saved_dir)
grn_inference.infer_grn()
```

### More
More info about LineageGRN can be found on our website https://lineagegrn.readthedocs.io. There you can find an API reference that give an introduction on how to use LineageGRN most effectively.
***

### Citation
If you find LineageGRN helpful for your research, please consider citing

