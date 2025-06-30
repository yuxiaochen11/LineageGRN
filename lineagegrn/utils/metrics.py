import numpy as np
import itertools
from sklearn.metrics import precision_recall_curve,roc_curve,auc
from .constant import * 
from lineagegrn.cell_fate_map import *

def compute_auroc(grn_infer, grn_true, threshold):
    """
    Evaluate the Area Under the Receiver Operating Characteristic (AUROC) for inferred and true GRNs.

    Parameters:
    grn_infer (DataFrame): A DataFrame containing inferred gene regulatory interactions.
    grn_true (DataFrame): A DataFrame containing true gene regulatory interactions.
    threshold (float): A threshold for determining the presence of an interaction.

    Returns:
    float: The calculated AUROC value.
    """
    # Flatten the DataFrames to 1D arrays
    grn_infer = grn_infer.values.flatten()
    grn_true = grn_true.values.flatten()
    
    # Convert inferred interactions to binary values based on the threshold
    for i in range(len(grn_infer)):
        if abs(grn_infer[i]) > threshold:
            grn_infer[i] = 1  # Interaction present
        else:
            grn_infer[i] = 0  # No interaction

    # Convert true interactions to binary values based on the threshold
    for i in range(len(grn_true)):
        if abs(grn_true[i]) > threshold:
            grn_true[i] = 1  # Interaction present
        else:
            grn_true[i] = 0  # No interaction

    # Calculate false positive rate and true positive rate
    fpr, tpr, thresholds = roc_curve(grn_true, grn_infer)
    
    # Calculate the area under the ROC curve
    auroc = auc(fpr, tpr)
    
    return auroc


def compute_auprc(grn_infer, grn_true, threshold):
    """
    Evaluate the Area Under the Precision-Recall Curve (AUPRC) for inferred and true GRNs.

    Parameters:
    grn_infer (DataFrame): A DataFrame containing inferred gene regulatory interactions.
    grn_true (DataFrame): A DataFrame containing true gene regulatory interactions.
    threshold (float): A threshold for determining the presence of an interaction.

    Returns:
    float: The calculated AUPRC value.
    """
    # Flatten the DataFrames to 1D arrays
    grn_infer = grn_infer.values.flatten()
    grn_true = grn_true.values.flatten()
    
    # Convert inferred interactions to binary values based on the threshold
    for i in range(len(grn_infer)):
        if abs(grn_infer[i]) > threshold:
            grn_infer[i] = 1  # Interaction present
        else:
            grn_infer[i] = 0  # No interaction

    # Convert true interactions to binary values based on the threshold
    for i in range(len(grn_true)):
        if abs(grn_true[i]) > threshold:
            grn_true[i] = 1  # Interaction present
        else:
            grn_true[i] = 0  # No interaction

    # Calculate precision and recall
    precision, recall, thresholds = precision_recall_curve(grn_true, grn_infer)
    
    # Calculate the area under the precision-recall curve
    auprc = auc(recall, precision)
    
    return auprc


def normalize_to_zero_mean(df):
    overall_mean = df.values.mean()  
    normalized_df = df - overall_mean 
    return normalized_df


def sort_dataframe_by_first_column(df):
    sorted_df = df.sort_values(by=df.columns[0], ascending=False)
    return sorted_df

def count_greater_than_threshold(df, threshold):
    count_series = ((df > threshold).sum(axis=1))/df.shape[1]
    return count_series.to_frame(name='count')

def remove_nan_elements(matrix1, matrix2):
    filtered_value= ~np.isnan(matrix1)
    filtered_matrix1 = matrix1[filtered_value]
    filtered_matrix2 = matrix2[filtered_value]
    filtered_matrix1[filtered_matrix1 < 0.1] = 0  
    filtered_matrix1[filtered_matrix1 > 0.1] = 1
    filtered_matrix2[filtered_matrix2 <0.1] = 0
    filtered_matrix2[filtered_matrix2 > 0.1] = 1
    return filtered_matrix1, filtered_matrix2

def match_ExpData(data_inf,data_real,rank):
    distance_matrix = np.zeros((len(data_inf), len(data_real))) 
    for i in range(len(data_inf)):
        for j in range(len(data_real)):
            filtered_data_inf=data_inf[i]
            filtered_data_real=data_real[j]
            filtered_data_inf,filtered_data_real=remove_nan_elements(data_inf[i], data_real[j])
            distance_matrix[i, j] =np.sum(filtered_data_inf != filtered_data_real)
    sums_list = []
    perm_list = []
    print(distance_matrix)
    for perm in itertools.permutations(range(len(data_inf))):
        total = sum(distance_matrix[i, perm[i]] for i in range(len(data_inf)))
        sums_list.append(total)
        perm_list.append(perm)
    sums_arr = np.array(sums_list)
    sorted_indices = np.argsort(sums_arr)
    sorted_sums = sums_arr[sorted_indices]
    min_val = sorted_sums.min()
    max_val = sorted_sums.max()
    norm_sums = (sorted_sums - min_val) / (max_val - min_val)
    idx = sorted_indices[rank - 1]
    loss_value = sums_arr[idx]
    norm_value = norm_sums[rank - 1]
    perm_choice = perm_list[idx]
    return norm_sums,loss_value, norm_value, perm_choice


def metric_perturbation_impact(grn_before_dict, grn_after_dict, threshold=0.05, top_k=30):
    """
    Compute average Jaccard index and Regulator overlap from two GRN dictionaries.

    Parameters:
    grn_before_dict: dict of DataFrame, keys are node names, values are GRNs (rows: targets, columns: regulators)
    grn_after_dict: dict of DataFrame, same structure as grn_before_dict
    threshold: float, threshold to binarize edges for Jaccard index
    top_k: int, number of top regulators to consider for overlap

    Returns:
    avg_jaccard: float, average Jaccard index across nodes
    avg_overlap: float, average regulator overlap across nodes
    """
    jaccard_scores = []
    overlap_scores = []

    common_keys = grn_before_dict.keys() & grn_after_dict.keys()

    for key in common_keys:
        grn_before = grn_before_dict[key].copy()
        grn_after = grn_after_dict[key].copy()

        # Ensure same shape and order
        grn_before = grn_before.loc[grn_after.index, grn_after.columns]

        # Jaccard index
        bin_before = grn_before.abs() > threshold
        bin_after = grn_after.abs() > threshold
        intersection = (bin_before & bin_after).values.sum()
        union = (bin_before | bin_after).values.sum()
        jaccard = intersection / union if union > 0 else 1.0
        jaccard_scores.append(jaccard)

        # Regulator overlap (columns are regulators)
        reg_strength_before = grn_before.abs().sum(axis=0)
        reg_strength_after = grn_after.abs().sum(axis=0)
        top_before = set(reg_strength_before.nlargest(top_k).index)
        top_after = set(reg_strength_after.nlargest(top_k).index)
        overlap = len(top_before & top_after) / top_k
        overlap_scores.append(overlap)

    avg_jaccard = sum(jaccard_scores) / len(jaccard_scores) if jaccard_scores else 0
    avg_overlap = sum(overlap_scores) / len(overlap_scores) if overlap_scores else 0

    return avg_jaccard, avg_overlap  

def match_ExpData(data_inf, data_real, rank=1, sample_size=10000):
    N = len(data_inf)
    D = np.zeros((N, N), dtype=float)
    for i in range(N):
        for j in range(N):
            inf, real = remove_nan_elements(data_inf[i], data_real[j])
            D[i, j] = np.sum(inf != real)

    if N <= 10:
        perms = list(itertools.permutations(range(N)))
    else:
        perms_set = {tuple(range(N))}
        while len(perms_set) < sample_size:
            perms_set.add(tuple(np.random.permutation(N)))
        perms = list(perms_set)

    losses = np.array([sum(D[i, p[i]] for i in range(N)) for p in perms])
    order = np.argsort(losses)
    sorted_losses = losses[order]
    sorted_perms = [perms[i] for i in order]

    mn, mx = sorted_losses[0], sorted_losses[-1]
    norm_sums = (sorted_losses - mn) / (mx - mn) if mx > mn else np.zeros_like(sorted_losses)

    k = min(rank, len(norm_sums)) - 1
    loss_value = sorted_losses[k]
    norm_value = norm_sums[k]
    perm_choice = sorted_perms[k]

    return norm_sums, loss_value, norm_value, perm_choice, sorted_perms