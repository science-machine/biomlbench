# BioML-bench

**THIS WAS GENERATED IN COLLABORATION WITH CLAUDE. EXPECT BUGS.**

See the docs here: [docs/index.md](docs/index.md)


## biomlbench core TODOs (unordered)

- [ ] Decide whether to enable multi-target problems in Polaris in the future (currently we just take the first target if more than one). 
- [ ] Add `experiments/splits` which will contain files like `medical-imaging.txt` that contain a list of task ids for only a given domain, or type of benchmarking run (e.g., lite, all).
- [ ] Decide if we want to have `experiments/familiarity` like the original MLE-bench. This would involve measuring the relationship between model performance and the "familiarity" of the model for the task description. This tries to assess if the model is "cheating" by already understanding the task.
- [ ] Need to decide if we can add human baselines in for Polaris. We would probably have to manually scrape the leaderboard (possibly with playwright) and then we would need to add a `supports_human_baselines` method to the `PolarisDataSource` class which returns `True`. And we would need to add a `get_human_baselines` method to the `PolarisDataSource` class which returns a pandas dataframe with the human baselines.
- [ ] Need an automated pipeline for scraping all the polaris and kaggle tasks and creating tasks/ directories for them.
- [ ] Revisit the `generate_submission` method. It currently assumes that the targets in train.csv are the last column. It also requires data in CSV format...
- [ ] Update the models used by the AIDE agent.
- [ ] Find a way to provide task-specific Docker containers for agents. Additionally, need to figure out whether we can enable more flexible base environment container (not hardcoding pytorch version, etc.)
- [ ] Figure out what they did in [the MLE-bench version of aideml](https://github.com/WecoAI/aideml/compare/main...thesofakillers:aideml:main) and see if we should implemnt that.
- [ ] Need to figure out what to do when no human baselines are available. How would scoring work? Also can we realistically combine scores across databases of benchmarks?
- [ ] Need to decide whether we want to incude obfuscated instructions/descriptions for tasks. This would probably take a lot of work but would be useful for benchmarking.
- [ ] Need to implement filtering for tasks by domain, type, complexity (e.g., only ADMET tasks, only imaging tasks, only lite tasks)
- [ ] Improve documentation of built-in agents.
- [ ] Decide if we're going to add in unit tests.
- [ ] Decide if we want to bother scraping subsequent polarihub leaderboard pages (not a lot of examples with more than one page).
- [ ] Decide if we want to have the option to let models see the public leaderboard (would be a more fair comparison to human performance).

Tasks:

- [ ] For `adaptyv-bio-egfr-binders-binary-cls-v0` we need to get the leaderboard from [this page](https://www.adaptyv.com/leaderboard/adaptyv-bio-egfr-binders-binary-cls-v0).


## How to wrap a new benchmark database

Example: Polaris.

To wrap polaris, we needed to:

1. Understand the data source. In this case, Polaris has a high-level API that can be used to automatically download and split data.
2. Create a new data source module in `biomlbench/data_sources/`. In the case of Polaris, we created `biomlbench/data_sources/polaris.py`.
3. Implement the `XYZDataSource` class. This class must implement the `download` and `get_leaderboard` methods.

