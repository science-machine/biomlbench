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
