# RxRx-compound-gene-activity-benchmark

![Recursion logo](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSgn5ii7g2UEVORKDmzzR0BKwzGr8Arz2ecTg&s)

# Recapitulating known compound-gene relationships 

This benchmark evaluates the **zero-shot prediction** of compound-gene activity.

Performance is measured separately for each compound using Average Precision and AUC-ROC. Per-compound measures are then aggregated into one global performance metric. 

 ## Maps of Biology and Chemistry
At Recursion, we build maps of biology and chemistry to explore uncharted areas of disease biology, unravel its complexity, and industrialize drug discovery. We use deep learning models to embed high dimensional representations of biology (e.g. phenomics, transcriptomics). This allows us to create representations that can be compared and contrasted to predict trillions of relationships across biology and chemistry — even without physically testing all of the possible combinations. Just as a map helps to navigate the physical world, our maps are designed to help us understand as much as we can about the connectedness of human biology so we can navigate the path to new medicines more efficiently.

We evaluated our recently released [`OpenPhenom-S/16`](https://www.rxrx.ai/phenom) model, as well as our proprietary `Phenom-1` and `Phenom-2` models. **Can you do better?**

![Benchmarking results of Phenom, OpenPhenom and the baselines](https://cdn.prod.website-files.com/5cb63fe47eb5472014c3dae6/6733a2d8e1ab553f6f26861e_compound_auc_roc_vs_mAP%20(1).png)

## Resources
This benchmark was released alongside [`RxRx3-Core`](https://polarishub.io/datasets/recursion/rxrx3-core) and [`OpenPhenom-S/16`](https://www.rxrx.ai/phenom).

 ### RxRx3-Core
A challenge dataset in phenomics optimized for the research community. RxRx3-core includes labeled images of 735 genetic knockouts and 1,674 small-molecule perturbations drawn from the RxRx3 dataset, image embeddings computed with OpenPhenom-S/16, and associations between the included small molecules and genes. The dataset contains 6-channel Cell Painting images and associated embeddings from 222,601 wells but is less than 18Gb, making it incredibly accessible to the research community.

### OpenPhenom-S/16
OpenPhenom-S/16is a foundation model that flexibly processes microscopy images into general-purpose embeddings. In other words, OpenPhenom-S/16 can take a series of microscopy channels and create a meaningful vector representation of the input image. This enables robust comparison of images, and other data science techniques to decode any biology or chemistry within such images. **The precomputed embeddings of OpenPhenom-S/16 are included in the RxRx3-core dataset!**

---

**Source:** [Polaris Hub - recursion/rxrx-compound-gene-activity-benchmark](https://polarishub.io)  
**Task Type:** single_task  
**Main Metric:** pr_auc{"aggregation": "mean", "default": null, "group_by": "treatment", "on_error": "ignore"}

## Data Format

This task uses the Polaris data source system:
- Data is downloaded from Polaris Hub using the PolarisDataSource
- Molecule column: `treatment`
- Target column: `{'active'}` (first target if multiple available)

## Files

- `train.csv`: Training data with molecules and targets
- `test_features.csv`: Test features with ID column
- `sample_submission.csv`: Example submission format

## Evaluation

Uses official Polaris evaluation system with benchmark `recursion/rxrx-compound-gene-activity-benchmark`.
Main metric: **pr_auc{"aggregation": "mean", "default": null, "group_by": "treatment", "on_error": "ignore"}**

## Source

Auto-generated from [Polaris Hub](https://polarishub.io/).
