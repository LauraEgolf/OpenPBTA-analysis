## Counting git contributions to a branch

### Purpose

The purpose of this module is to count the git contributions to this repo in two ways:

1. Count the number of analysis modules – directories in `analyses/` that are included in the manuscript, specifically – that individual contributors have authored commits in. 
  You can find that information in `results/module_contribution_counts.tsv` and a _less summarized_ version that lists all contributors to a module in `results/module_contributors.tsv`.
2. Count the number of commits to the entire repository individual contributors have authored.

### Steps

The shell script `run-count-contributions.sh` runs the following steps:

* `01-count-contributions.sh` which is a shell script that generates intermediate TXT files (in `scratch/count-contributions`) with the output of `git shortlog` and `git log`.
* `02-format-contributions.Rmd` which turns the output of `git shortlog` and `git log` into readable files. 
This notebook additionally collapses people that have their contributions divided across multiple names by the logs and is responsible for removing the analysis modules that are not included in the manuscript. 
See the notebook for more details!

### GitHub Actions workflow

To keep these stats reasonably up-to-date, we use a GitHub Actions (GHA) workflow ([`.github/workflows/count-git-contributions.yml`](https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/.github/workflows/count-git-contributions.yml)) to rerun the module every Wednesday at 14:00 UTC.

The workflow runs `bash run-count-contributions.sh` on the `master` branch and files a new pull request with the updated tables and HTML file from the notebook.
This PR is intended to be reviewed by an organizer for accuracy before merging.

The GHA workflow is the main way we expect the module to be used.

### Running the module manually

To count git contributions on :warning: your current branch :warning:, you can run the following:

```sh
bash run-count-contributions.sh
```

#### Checking out the correct branch

**You most likely want to count contributions to the `AlexsLemonade` `master` branch!**
We include some instructions for doing this locally below.

`AlexsLemonade` should be set as an `upstream` remote, which you can check with:

```sh
git remote -v
```

Which should include the following if `AlexsLemonade` is set as `upstream`:

```sh
upstream	https://github.com/AlexsLemonade/OpenPBTA-analysis.git (fetch)
upstream	https://github.com/AlexsLemonade/OpenPBTA-analysis.git (push)
```

You likely have your own `master` branch locally, so to checkout `upstream/master` with a new local name (using `upstream-master` below), you can use the following:

```sh
git checkout -b upstream-master upstream/master
```

If you already have a local `upstream-master` branch, check it out and make sure it's up to date!

Now you're ready to run the module!

#### Committing the changes from a manual run

Once you've run `bash run-count-contributions.sh` to generate stats for the `AlexsLemonade/master` branch, you'll need to _create a new branch_ and commit the modified files to that new branch.

