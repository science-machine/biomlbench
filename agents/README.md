# Agents

**See documentation on [creating agents](../docs/developer/creating_agents.md)**

MLE-Bench original documentation:

We evaluated three open-source agents. We slightly modified each agent to elicit better capabilities. We also developed a "dummy" agent, used to check that the environment is configured correctly. Each agent, alongside the link to our fork, is listed below. Each agent has an associated ID which we use to identify it within our repo.

Agent | ID | Fork
:-----:|:-----:|:-----:
dummy | dummy | N/A
[AIDE](https://www.weco.ai/blog/technical-report) | aide | https://github.com/AIxBio/aideml
[MLAgentBench](https://openreview.net/forum?id=1Fs1LvjYQW) | mlagentbench | https://github.com/JunShern/MLAgentBench
[OpenHands](https://arxiv.org/abs/2407.16741) | opendevin | https://github.com/thesofakillers/OpenHands

## Prerequisites
If you want to run these agents locally:
- Install [Docker](https://docs.docker.com/engine/install/)
- Install [Sysbox](https://github.com/nestybox/sysbox). See [Security](#Security) below for more information
- (Optional) Install [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) to run agents with GPUs

## Getting Agent Images

### Option 1: Use Prebuilt Images (Recommended)

Pull prebuilt images and tag them for local use:

```bash
# Pull all prebuilt images and tag them locally
./scripts/pull_prebuilt_images.sh
```

This is much faster than building locally and pulls the following:
- `millerh1/biomlbench-env:v0.1a` → tagged as `biomlbench-env`
- `millerh1/aide:v0.1a` → tagged as `aide`
- `millerh1/biomni:v0.1a` → tagged as `biomni`
- `millerh1/mlagentbench:v0.1a` → tagged as `mlagentbench`
- `millerh1/stella:v0.1a` → tagged as `stella`
- `millerh1/dummy:v0.1a` → tagged as `dummy`

### Option 2: Build Images Locally

To build an image for an agent with ID `<agent>`, run:

```bash
# Build base environment first
./scripts/build_base_env.sh

# Then build specific agent
./scripts/build_agent.sh <agent>

# Or use the manual build command:
export SUBMISSION_DIR=/home/submission
export LOGS_DIR=/home/logs
export CODE_DIR=/home/code
export AGENT_DIR=/home/agent

docker build --platform=linux/amd64 -t <agent> agents/<agent>/ --build-arg SUBMISSION_DIR=$SUBMISSION_DIR --build-arg LOGS_DIR=$LOGS_DIR --build-arg CODE_DIR=$CODE_DIR --build-arg AGENT_DIR=$AGENT_DIR
```

## Running agents

The `biomlbench run-agent` command allows you to run agents locally on a given set of tasks. In the `experiments/splits/` directory, we have several files, each containing a set of task IDs. The `experiments/splits/all.txt` file contains all tasks. The `experiments/splits/caco2-wang.txt` split just contains the caco2-wang task, which is useful for testing. For example, to run the dummy agent on the caco2-wang task, you can run:

```console
biomlbench run-agent --agent dummy --task-list experiments/splits/caco2-wang.txt
```

Running `biomlbench run-agent` creates a "run group" directory in the `runs/` directory. The run group directory will contain a subdirectory for each task that the agent was evaluated on, containing the agent's logs, code, and submission. A `metadata.json` file and `submission.jsonl` file are automatically created when the run finishes. You can then grade this run directly:

```bash
biomlbench grade --submission runs/<run-group>/submission.jsonl --output-dir runs/<run-group>
```

If you'd like to update the configuration of the container, you can edit the default container config in `environment/config/container_configs/default.json`, or specify a custom container config JSON file when executing `biomlbench run-agent` with the `--container-config` flag. If you'd like to run the agent with a GPU, you can set `"gpus": -1` in the container config JSON file.

### dummy agent
We used the dummy agent to verify that our environment was configured correctly. It performs the following checks:

1. Prints to stdout whether its python script `main.py` is running with root access
2. Prints to stdout the Python interpreter the dummy agent is using
3. Attempts to use the sample submission of the current competition as it's submission. Note, this sometimes fails as sample submissions are sometimes not in the expected format. For example, some competitions compress the sample submission along with the other data into a single file
4. Checks that it can't read the "private" data; this includes the labels of the test set
5. Checks it has read/write access to `/home`, it's working directory

### Security

By default we run AIDE, MLAB, and the dummy agent using the [sysbox](https://github.com/nestybox/sysbox) runtime. We recommend users to use this runtime, or other runtimes with sufficient security guarantees, when running models on this benchmark. Since container runtimes like Docker can share host resources with containers, a malicious container could potentially compromise a host; sysbox is designed to mitigate this risk by providing enhanced isolation between containers and the host system

OpenHands uses it's own docker containers during execution, resulting in docker-in-docker (DinD). The sysbox runtime does not yet support GPU-passthrough in DinD ([link](https://github.com/nestybox/sysbox/issues/50)). We were therefore required to run OpenHands in "privileged" mode, which we deemed acceptable given that our (internal) evaluation infrastructure does not rely on the Docker container-to-host isolation or the virtual machine boundary for security. We wish users to acknowledge this security risk, thus require users to set the `I_ACCEPT_RUNNING_PRIVILEGED_CONTAINERS` environment variable to `True` when running OpenHands through `biomlbench run-agent`.
