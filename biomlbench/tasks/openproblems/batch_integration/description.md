# batch_integration

## Summary
Remove unwanted batch effects from scRNA data while retaining biologically meaningful variation.

## Description
In this task we evaluate batch integration methods on their ability to remove batch effects in the data while conserving variation attributed to biological effects.
As input, methods require either normalised or unnormalised data with multiple batches and consistent cell type labels.
The batch integrated output can be a feature matrix, a low dimensional embedding and/or a neighbourhood graph.
The respective batch-integrated representation is then evaluated using sets of metrics that capture how well batch effects are removed and whether biological variance is conserved.
We have based this particular task on the latest, and most extensive benchmark of single-cell data integration methods.


## Motivation
As single-cell technologies advance, single-cell datasets are growing both in size and complexity.
Especially in consortia such as the Human Cell Atlas, individual studies combine data from multiple labs, each sequencing multiple individuals possibly with different technologies.
This gives rise to complex batch effects in the data that must be computationally removed to perform a joint analysis.
These batch integration methods must remove the batch effect while not removing relevant biological information.
Currently, over 200 tools exist that aim to remove batch effects scRNA-seq datasets [@zappia2018exploring].
These methods balance the removal of batch effects with the conservation of nuanced biological information in different ways.
This abundance of tools has complicated batch integration method choice, leading to several benchmarks on this topic [@luecken2020benchmarking; @tran2020benchmark; @chazarragil2021flexible; @mereu2020benchmarking].
Yet, benchmarks use different metrics, method implementations and datasets. Here we build a living benchmarking task for batch integration methods with the vision of improving the consistency of method evaluation.


---

**Source:** [OpenProblems - batch_integration](https://openproblems.bio)  
**Task Type:** benchmark_results  
**Data Type:** Results-only (pre-computed)

## Available Data

This task provides access to pre-computed benchmark results from OpenProblems:

- **Methods**: 4 integration methods
- **Datasets**: 10 single-cell datasets
- **Metrics**: 8 evaluation metrics
- **Total Results**: 183 method-dataset combinations

## Evaluation Metrics

OpenProblems uses multiple metrics to evaluate batch integration:
- **KBET**: k-Nearest Neighbor Batch Effect Test
- **ASW (batch)**: Average Silhouette Width for batch mixing
- **ASW (label)**: Average Silhouette Width for label separation  
- **Graph connectivity**: Preservation of local neighborhoods
- **ARI**: Adjusted Rand Index
- **NMI**: Normalized Mutual Information
- **Isolated label ASW**: Silhouette width for isolated labels
- **Isolated label F1**: F1 score for isolated label detection

## Usage

This is a **results-only** task that displays pre-computed benchmarks. For actual submission and evaluation, visit [OpenProblems](https://openproblems.bio).

## Source

Auto-generated from [OpenProblems](https://openproblems.bio/).
