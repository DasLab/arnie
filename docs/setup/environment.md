# Setting up an arnie environment

`arnie` is a Python package to simplify interacting with various RNA prediction and analysis libraries. To work, the `arnie` package needs to know the location of those libraries on the local filesystem. `arnie` uses environment variables to point to package locations. 

## Environment Variables
Here we assume you've already installed a package you want arnie to use (visit the [supported packages page](/setup/packages.md) for more details about specific package installation requirements). Arnie expects environment variables in the form of "{package_name}_PATH". So for `contrafold`, we specify its installed location for arnie with `export contrafold_PATH=/path/to/executable/`. Certain packages require additional resources for arnie to operate. For example, `SpotRNA` also requires a pointer to the conda environment it is installed with. The [supported packages page](/setup/packages.md) details each package's expected environment variables.

Arnie also expects an `arnie_TMP` environment variable to define where arnie should write temporary files to. Some predictor packages write to files to generate their output; arnie uses the `arnie_TMP` location to support these packages.

## Arnie File
As a fallback, you can also specify an "arnie_file.txt" that defines these paths. There is an example arnie_file.txt included in the arnie repo that demonstrates the expected syntax. If using the arnie_file approach, you need to set an `ARNIEFILE` environment variable pointing to your arnie_file.txt (e.g, `export ARNIEFILE="/path/to/arnie/<my_file.txt>"`)

## Conda Environments
We recommend using [conda](https://anaconda.org/anaconda/conda) to set up private Python execution environments for your arnie operations. Conda simplifies the sometimes complicated process of managing Python dependencies by creating virtual environments that isolate installed packages. Conda also supports simplified distribution of a wide range of scientific Python libraries, and even a number of RNA structure packages. We recommend the following setup for your RNA science conda environment.
```
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
[Bioconda](https://bioconda.github.io/) and [Conda-Forge](https://conda-forge.org/) are distribution channels for conda packages. Bioconda, for instance, hosts `ViennaRNA` and `Eternafold` RNA packages

We set up an example conda environment to support our arnie work below. First we create the environment. Next we activate the environment, which sets up our isolated Python execution environment. After activation, we pip install arnie (which will be installed into the isolated environment with proper PYTHONPATH handling). 
```
conda create -n rna-env
conda activate rna-env
pip install arnie
```