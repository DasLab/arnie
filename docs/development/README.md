# Contributing

## Installing via Github
The `arnie` package source code is hosted on [Github](https://github.com/DasLab/arnie). You can clone the repo as below.

```
git clone https://github.com/DasLab/arnie.git
```

You can also use pip to install arnie from our Github repo:
```
pip install git+https://github.com/DasLab/arnie
```
This is particularly useful for testing new features internally before releasing on PyPI.

## Repo Organization

`src/arnie`: source code for the arnie package.

`docs`: docsify-based markdown documentation for the arnie package.

`tests`: unit tests 

`notebooks`: example jupyter notebooks with usage.

`scripts`: scripts for processing sequences in batch.

`parameter_files`: dir of various parameter files for packages, put here out of convenience.



## Github Issues
We use [Github issues](https://github.com/DasLab/arnie/issues) to coordinate development tasks and track feature development and bug fixes. If you run into problems while using `arnie`, please file an issue so that we can address the bug. Similarly, if you have a feature idea that could simplify your research, file an issue detailing your proposed feature. 

## Package Testing
Tests are located in the `tests` directory of the repo. We use the [pytest](https://docs.pytest.org/en/stable/) testing framework. Tests are run in the repo root directory. 

To run all the tests,
```
pytest
```
To run a specific test,
```
pytest tests/test_structure_handling.py
```
If you add new features or fix a bug, make sure to update the tests appropriately. 

## Package Distribution
We distribute arnie via the [Python Package Index](https://pypi.org/). The DasLab has a [PyPI account](https://pypi.org/user/daslab/) for all our packages, with `arnie` available [here](https://pypi.org/project/arnie/)

Arnie package release is automated via Github Actions. The [release workflow](https://github.com/DasLab/arnie/actions/workflows/release.yml) builds the package for distribution, publishes to PyPI and releases a Github release. The action is triggered on new git tag push. 

To push a new release, update the `pyproject.toml` version number as appropriate (we follow the [semantic versioning](https://semver.org/) standard). Next, define a matching git tag for the version number, and then push to Github.
```
git checkout master
git tag -a v1.1.0 -m "Arnie Release v1.1.0"
git push origin tag v1.1.0 
```