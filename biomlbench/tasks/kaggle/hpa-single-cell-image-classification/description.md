# Human Protein Atlas - Single Cell Classification

## Description

There are billions of humans on this earth, and each of us is made up of trillions of cells. Just like every individual is unique, even genetically identical twins, scientists observe differences between the genetically identical cells in our bodies.

Differences in the location of proteins can give rise to such cellular heterogeneity. Proteins play essential roles in virtually all cellular processes. Often, many different proteins come together at a specific location to perform a task, and the exact outcome of this task depends on which proteins are present. As you can imagine, different subcellular distributions of one protein can give rise to great functional heterogeneity between cells. Finding such differences, and figuring out how and why they occur, is important for understanding how cells function, how diseases develop, and ultimately how to develop better treatments for those diseases.

The study of a single cell enables the discovery of mechanisms too difficult to see with multi-cell research. The importance of studying single cells is reflected in the ongoing revolution in biology centered around technologies for single cell analysis. Microscopy offers an opportunity to study differences in protein localizations within a population of cells. Current machine learning models for classifying protein localization patterns in microscope images gives a summary of the entire population of cells. However, the single-cell revolution in biology demands models that can precisely classify patterns in each individual cell in the image.

The Human Protein Atlas is an initiative based in Sweden that is aimed at mapping proteins in all human cells, tissues, and organs. The data in the Human Protein Atlas database is freely accessible to scientists all around the world that allows them to explore the cellular makeup of the human body. Solving the single-cell image classification challenge will help us characterize single-cell heterogeneity in our large collection of images by generating more accurate annotations of the subcellular localizations for thousands of human proteins in individual cells.

This is a **weakly supervised multi-label classification problem**. Given images of cells from microscopes and labels of protein location assigned together for all cells in the image, you will develop models capable of segmenting and classifying each individual cell with precise labels.

## Task Overview

You are predicting **protein organelle localization labels for each cell** in the image. Border cells are included when there is enough information to decide on the labels.

### The Challenge: Weak Supervision

The labels you will get for training are **image level labels** while the task is to predict **cell level labels**. That is to say, each training image contains a number of cells that have collectively been labeled as described above and the prediction task is to look at images of the same type and predict the labels of each individual cell within those images.

As the training labels are a collective label for all the cells in an image, it means that each labeled pattern can be seen in the image but not necessarily that each cell within the image expresses the pattern. This imprecise labeling is what we refer to as **weak**.

During the challenge you will both need to **segment the cells** in the images and **predict the labels** of those segmented cells.

## Data Format

### Image Structure
Each sample consists of **four files** representing different filters on the subcellular protein patterns:

- **Red**: Microtubule channels  
- **Blue**: Nuclei channels
- **Yellow**: Endoplasmic Reticulum (ER) channels
- **Green**: Protein of interest (target for prediction)

The format is `[filename]_[filter color].png` for the PNG files.

### Image Specifications
- **Training images**: Available in .tif format
- **Test images**: Available in .png format  
- **Sizes**: Mix of 1728×1728, 2048×2048, and 3072×3072 pixels
- **Imaging modality**: Confocal microscopy (highly standardized)
- **Cell types**: 17 different cell types with highly different morphology

### Subcellular Localization Labels

There are **19 different labels** present in the dataset (18 labels for specific locations, and label 18 for negative and unspecific signal):

0. Nucleoplasm  
1. Nuclear membrane   
2. Nucleoli   
3. Nucleoli fibrillar center   
4. Nuclear speckles   
5. Nuclear bodies   
6. Endoplasmic reticulum   
7. Golgi apparatus   
8. Intermediate filaments  
9. Actin filaments  
10. Microtubules      
11. Mitotic spindle   
12. Centrosome   
13. Plasma membrane   
14. Mitochondria   
15. Aggresome   
16. Cytosol   
17. Vesicles and punctate cytosolic patterns   
18. Negative  

## Evaluation

Submissions are evaluated by computing **mAP (mean Average Precision)**, with the mean taken over the 19 segmentable classes of the challenge. It is identical to the OpenImages Instance Segmentation Challenge evaluation.

**Segmentation is calculated using IoU with a threshold of 0.6.**

## Submission File

For each image in the test set, you must predict a list of instance segmentation masks and their associated detection score (Confidence). The submission CSV file uses the following format:

```
ImageID,ImageWidth,ImageHeight,PredictionString
ImageAID,ImageAWidth,ImageAHeight,LabelA1 ConfidenceA1 EncodedMaskA1 LabelA2 ConfidenceA2 EncodedMaskA2 ...
ImageBID,ImageBWidth,ImageBHeight,LabelB1 ConfidenceB1 EncodedMaskB1 LabelB2 ConfidenceB2 EncodedMaskB2 …
```

**Note**: A mask MAY have more than one class. If that is the case, predict separate detections for each class using the same mask.

### Sample Submission
```
ID,ImageWidth,ImageHeight,PredictionString
721568e01a744247,1118,1600,0 0.637833 eNqLi8xJM7BOTjS08DT2NfI38DfyM/Q3NMAJgJJ+RkBs7JecF5tnAADw+Q9I
7b018c5e3a20daba,1600,1066,16 0.85117 eNqLiYrLN7DNCjDMMIj0N/Iz9DcwBEIDfyN/QyA2AAsBRfxMPcKTA1MMADVADIo=
```

### Mask Encoding

The binary segmentation masks are:
1. **RLE encoded** using COCO's mask API  
2. **zlib compressed** (RFC1950)
3. **base64 encoded** to be used in text format as EncodedMask

Example encoding function:
```python
import base64
import numpy as np
from pycocotools import _mask as coco_mask
import typing as t
import zlib

def encode_binary_mask(mask: np.ndarray) -> t.Text:
    """Converts a binary mask into OID challenge encoding ascii text."""
    
    # check input mask --
    if mask.dtype != np.bool:
        raise ValueError(
            "encode_binary_mask expects a binary mask, received dtype == %s" %
            mask.dtype)
    
    mask = np.squeeze(mask)
    if len(mask.shape) != 2:
        raise ValueError(
            "encode_binary_mask expects a 2d mask, received shape == %s" %
            mask.shape)
    
    # convert input mask to expected COCO API input --
    mask_to_encode = mask.reshape(mask.shape[0], mask.shape[1], 1)
    mask_to_encode = mask_to_encode.astype(np.uint8)
    mask_to_encode = np.asfortranarray(mask_to_encode)
    
    # RLE encode mask --
    encoded_mask = coco_mask.encode(mask_to_encode)[0]["counts"]
    
    # compress and base64 encoding --
    binary_str = zlib.compress(encoded_mask, zlib.Z_BEST_COMPRESSION)
    base64_str = base64.b64encode(binary_str)
    return base64_str
```

## Timeline

- **January 26, 2021**: Start Date
- **May 4, 2021**: Entry deadline
- **May 4, 2021**: Team Merger deadline
- **May 11, 2021**: Final submission deadline

## Prizes

- **1st Place**: $12,000
- **2nd Place**: $8,000  
- **3rd Place**: $5,000

## Citation

Casper Winsnes, Emma Lundberg, Maggie, Phil Culliton, Trang Le, UAxelsson, and Wei Ouyang. Human Protein Atlas - Single Cell Classification. https://kaggle.com/competitions/hpa-single-cell-image-classification, 2021. Kaggle.

# Data

## Dataset Description

The goal of this competition is to detect and classify protein organelle localization labels for each cell in fluorescent microscopy images. The challenge involves both cell segmentation and multi-label classification with weak supervision.

### What is meant by weak image-level labels?

The labels you will get for training are image level labels while the task is to predict cell level labels. That is to say, each training image contains a number of cells that have collectively been labeled as described above and the prediction task is to look at images of the same type and predict the labels of each individual cell within those images.

As the training labels are a collective label for all the cells in an image, it means that each labeled pattern can be seen in the image but not necessarily that each cell within the image expresses the pattern. This imprecise labeling is what we refer to as weak.

During the challenge you will both need to segment the cells in the images and predict the labels of those segmented cells.

### Files

- **train/** - training images (in .tif format)
- **test/** - test images (in .png format)  
- **train.csv** - filenames and image level labels for training set
- **sample_submission.csv** - filenames for test set and submission format guide

### Data Columns

- **ID**: Base filename of the sample (each sample has four color files: red, blue, yellow, green)
- **Label**: Image level labels (training) / cell level labels (submission)

### Dataset Size

- **Total files**: 89,542 files
- **Dataset size**: 170.03 GB
- **File types**: png, tfrec, csv

# Required Output

- `submission.csv`: Submission file. Must be in the same format as `sample_submission.csv`. 