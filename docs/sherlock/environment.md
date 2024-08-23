# Sherlock

[Sherlock](https://www.sherlock.stanford.edu/), Stanford's high-performance computing cluster, is a useful resource for compute-intensive arnie tasks. 

If you've never worked with cluster computing before, there are some differences in how things work from your laptop. The [Sherlock docs](https://www.sherlock.stanford.edu/docs/) are a great place to start, and the Sherlock team even does [onboarding sessions and office hours](https://www.sherlock.stanford.edu/docs/#onboarding-sessions) to help new users.

## Storage on Sherlock
The first thing to understand when setting up your Sherlock environment is where to store data. Sherlock offers several [data storage systems](https://www.sherlock.stanford.edu/docs/storage/) tailored for specific needs. We recommend using them as follows:
- `$HOME`: storage for your rna-environment, miniconda install, source code, etc
- `$GROUPHOME`: storage for shared resources or projects that other lab members may access
- `$SCRATCH/$GROUP_SCRATCH`: high-performance storage for large datasets and temporary files (WARNING: files on SCRATCH and GROUP_SCRATCH are automatically purged
 90 days after their last content modification; make sure you back up data there before this window)
- `$LSCRATCH`: node-local SSD; useful for specific jobs where high IOPS are important or when performing large batch jobs that may impact group resources
- `$OAK`: long-term storage of large research datasets

## Installing software on Sherlock
Sherlock provides specific scientific computing software pre-installed on the Sherlock system via ["modules"](https://www.sherlock.stanford.edu/docs/software/modules/). These modules are selected and maintained to provide maximum compatibility and reduce dependency conflicts. You can search for available modules [here](https://www.sherlock.stanford.edu/docs/software/list/). If you want to use Sherlock's module system instead of setting up your own Python environment, we recommend the following modules to load for arnie (and various downstream prediction algorithms):
```
module load python/3.6.1
module load py-numpy/1.18.1_py36 
module load py-pandas/1.0.3_py36
module load py-scipy/1.4.1_py36
module load gcc
module load glpk
module load mpfr
```

## Setting up your environment on Sherlock
Setting up your environment on Sherlock is fairly straightforward. First, set up a folder for yourself in $GROUPHOME (`mkdir $GROUPHOME/{your_name}`) to store shared resources. Next, we'll install [miniconda](https://docs.anaconda.com/miniconda/) for Python environment management and package installation. After installing miniconda, we can configure conda with useful package channels, create and activate an rna-env environment to store our packages, and install arnie. We provide an environment.yaml folder at `$GROUPHOME/rna-env/rna-environment.yaml` to create a standard environment with some standard packages. 
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n rna-env -f $GROUPHOME/rna-env/rna-environment.yaml
conda activate rna-env
pip install arnie
```

Your new conda environment has arnie and a few predictors installed. However, many prediction libraries are not available via conda or pip and usually require installing from source. The lab maintains a directory of predictors on Sherlock that you should copy to your $HOME directory.
```
cd $HOME
git clone $GROUPHOME/rna-env
```
Predictors are stored under `rna-env/predictors`. If you add new predictors in the course of your work, make sure to push your updates back to the $GROUPHOME origin repo.

Now that your environment is set up, let's take a look at [using Sherlock for compute jobs with arnie](jobs.md).