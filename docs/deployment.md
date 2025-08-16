# Google cloud deployment with biomlbench

We provide scripts for deploying agents on GCP. You must have a `google-application-credentials.json` file.

Add the following 3 fields to your `.env` file:

```bash
GOOGLE_APPLICATION_CREDENTIALS=/path/to/google-application-credentials.json
GCP_PROJECT_ID=...  # can be found in google-application-credentials.json
GCP_ZONE=...  # up to you!
```

## Google cloud CLI installation

Follow the instructions below to install the google cloud CLI, which is needed to run deployment.

### MacOS (Homebrew)

```bash
brew install --cask google-cloud-sdk
```

### Linux

```bash
curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-468.0.0-linux-x86_64.tar.gz
tar -xf google-cloud-cli-468.0.0-linux-x86_64.tar.gz

# Run the installer
./google-cloud-sdk/install.sh

# Restart your shell or source the updated profile
source ~/.bashrc   # or ~/.zshrc
```

### Windows

Easiest way is to download the [installer](https://cloud.google.com/sdk/docs/install#windows) directly from the google cloud SDK website.

## Setting up google cloud

```bash
gcloud auth activate-service-account --key-file=/path/to/google-application-credentials.json
gcloud config set project <PROJECT_ID>
```

Again the `PROJECT_ID` can be grabbed from the JSON file.

## Deploying biomlbench on a google cloud VM

For running biomlbench on google cloud you should use the pre-baked OS image that contains all the prebuilt agent docker images and datasets. To create an image with a standard L4 GPU:

```bash
gcloud compute instances create biomlbench --zone=us-central1-a --machine-type=g2-standard-8 --maintenance-policy=TERMINATE --image=biomlbench --boot-disk-size=500G
```

All of the ProteinGym, single-cell, and PolarisHub datasets have been downloaded and prepared, as well as the following four kaggle tasks:
    - `histopathologic-cancer-detection`
    - `osic-pulmonary-fibrosis-progression`
    - `rsna-miccai-brain-tumor-radiogenomic-classification`
    - `uw-madison-gi-tract-image-segmentation`

Where `--machine-type=g2-standard-8` specifies the machine type with an L4 GPU. For the tasks that require large-scale deep learning models you can specify `--machine-type=a2-highgpu-1g`, which gives an instance with a single 40Gb A100 GPU.

Then you can `ssh` into the VM as follos:

```bash
gcloud compute ssh runner@biomlbench
```

You can also run a specific command:

```bash
gcloud compute ssh runner@biomlbench --command="cd biomlbench && source .venv/bin/activate && biomlbench run-agent --agent aide --task-id proteingym-dms/A0A1I9GEU1_NEIME_Kennouche_2019"
```

All tasks from ProteinGym and Polaris have been prepared, as have the four kaggle tasks in the google doc.
