# Setting up Arnie (Mac and Linux)

Arnie is a set of wrappers to call other packages. So in order for Arnie to work right, we have to do a little bit of lifting at the start to 1) download and compile those packages, and 2) tell Arnie where they are in your arnie file.

You don't have to install all of the packages described here for Arnie to work. Your needs will dictate what other packages you set up.

- More accurate thermodynamic predictions? Set up EternaFold

- SHAPE-directed structure prediction? Set up Vienna or RNAstructure

- Benchmarking packages / skeptic of all?  Set up all of them!

This document contains instructions to download, compile where needed, and set paths for all the major packages and versions that Arnie supports. This document does not describe the python syntax -- for that, go to `notebooks/start_here.ipynb`.

### Point python to arnie location and define Arnie file location

1. Add Arnie location to your python path in your .bashrc, i.e.

```
export PYTHONPATH=$PYTHONPATH:/path/to/arnie
```

2. Create a variable in your .bashrc:

```
export ARNIEFILE="/path/to/arnie/my_arnie_file.txt"
```
Next we're going to create that arnie file.

### Creating the Arnie file

At the end of the day you'll have an example arnie file that looks like this:

```
# Vienna RNAfold 2 Mac using installer:
vienna_2: /usr/local/bin

# NUPACK build:
nupack: /path/to/nupack3.0.6/src/bin

# CONTRAfold build
contrafold_2: /path/to/contrafold-se/src

# RNASoft build
rnasoft: /path/to/MultiRNAFold

# RNAstructure (precompiled command line src)
rnastructure: /path/to/RNAstructure/exe

# Path to EternaFold
eternafold: /path/to/EternaFold/src

# (DEPRECATED, better to put path to Eternafold)
# Path to EternaFold parameters
eternafoldparams: /path/to/EternaFold/parameters/EternaFoldParams.v1

# LinearFold build
linearfold: /path/to/LinearFold/bin

#LinearPartition build
linearpartition: /path/to/LinearPartition/bin

#directory to write temp files
TMP: /tmp
```
## EternaFold (`package='eternafold'`)

Go to [https://eternagame.org/about/software](https://eternagame.org/about/software). Follow the instructions under "Request License" under EternaFold to download.

![](doc_ASSETS/Untitled%203.png)

Check the build by running:

`./src/contrafold predict test.seq --params parameters/EternaFoldParams.v1`

Output:
```
Training mode: 
Use constraints: 0
Use evidence: 0
Predicting using MEA estimator.
>test.seq
CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG
>structure
(((((((((((((......))))))..)....((((.....))))...))))))
```

Set in your arnie file:

```
#Path to EternaFold
eternafold: /path/to/EternaFold/src
```

EternaFold can also be set up to use the parameters with LinearFold/LinearPartition. For more information, see documentation in [EternaFold](https://eternagame.org/about/software).

## Vienna RNAFold

### Existing binaries:

(I use this for Mac.)

Download a precompiled binary from [https://www.tbi.univie.ac.at/RNA/#download](https://www.tbi.univie.ac.at/RNA/#download) .

![](doc_ASSETS/Untitled.png)

Click on the icon of your OS to download the installer, and click through the installer. This should install the Vienna executables to `/usr/local/bin`. Check by navigating there and trying to run `RNAfold -h`.

Then set in the arnie file:

```
# Vienna RNAfold 2 Mac using installer:
vienna_2: /usr/local/bin
```

### Compile from source:

(I use this for linux.)

Download the source code from [https://www.tbi.univie.ac.at/RNA/#download](https://www.tbi.univie.ac.at/RNA/#download) :

![](doc_ASSETS/Untitled%202.png)

```
gunzip ViennaRNA-2.4.14.tar.gz
tar -xvf ViennaRNA-2.4.14.tar
cd ViennaRNA-2.4.14
./configure --prefix /path/to/vienna/bin # make this where you want the executables to be installed
make install
```

I get errors about not being able to find `EXTERN.H` from perl wrappers, but the executable builds still work.

Check the build by navigating to `/path/to/vienna/bin` and running `RNAfold -h`.

Then set in the arnie file:

```
# Vienna RNAfold 2 Linux build from source:
vienna_2: /path/to/ViennaRNA-2.4.14/src/bin
```

## NUPACK (`package='nupack'`)

Download from [http://www.nupack.org/downloads](http://www.nupack.org/downloads). (Requires registering to sign a license agreement).

Notes: I have found that NUPACK 3.2.2 does not work on Mac, and instead use NUPACK 3.0.6 on Mac. Nupack 3.2.2 works on linux though. 

```
cd nupack3.0.6
make
```

Check build by running `./bin/mfe -h`.

Then set in arnie file:

```
# NUPACK build:
nupack: /path/to/nupack3.0.6/src/bin
```

## CONTRAfold 2 (`package='contrafold'`)

I recommend downloading the code that is at [https://github.com/csfoo/contrafold-se](https://github.com/csfoo/contrafold-se).  Even though this is named CONTRAfold-SE, the repo itself does not contain any of the parameter files from the CONTRAfold-SE project, and the authors cleaned up some build issues. If you run this code with its default parameters, you are running CONTRAfold 2.

```
git clone https://github.com/csfoo/contrafold-se.git
cd contrafold-se/src
make
```

If you get the following error:
```
Utilities.cpp:342:17: error: use of undeclared identifier 'mkdtemp'
```
Try adding `#include <unistd.h>` to `Utilities.hpp`.

Check the build by running `./src/contrafold -h`.

Set in your arnie file:

```
# CONTRAfold build
contrafold_2: /path/to/contrafold-se/src
```

If you get the error:

```
‘LONG_MIN’ was not declared in this scope
```

Try adding `#include <climits>` to Utilities.hpp.

## RNAStructure (`package='rnastructure'`)

Download RNAstructure (pre-compiled command-line interface) for your OS from [https://rna.urmc.rochester.edu/RNAstructure.html](https://rna.urmc.rochester.edu/RNAstructure.html).

![](doc_ASSETS/Untitled%204.png)

Unzip the tarball and test the executables in `RNAstructure/exe` by running `./partition -h`.

Add to Arnie file:

```
# RNAstructure (precompiled command line src)
rnastructure: /path/to/RNAstructure/exe
```

## RNASoft (`package='rnasoft'`)

Download the MultiRNAFold package from [http://www.rnasoft.ca/download.html](http://www.rnasoft.ca/download.html).  (Update 2022: Link no longer is functional.)

(Below tested for version 2.1, 2018).

Unzip the MultiRNAfold tarball.

For Mac:

May need to copy this `malloc.h` into the `include` dir: [https://opensource.apple.com/source/Libc/Libc-825.26/include/malloc/malloc.h.auto.html](https://opensource.apple.com/source/Libc/Libc-825.26/include/malloc/malloc.h.auto.html).

Then type

```
make depend
make
```

To check build, try `./simfold -h`. Then put in arnie file:

```
#RNASoft build
rnasoft: /path/to/MultiRNAFold
```

There are additional parameter files for other types of training that can be downloaded at [https://www.cs.ubc.ca/~andrones/param-training/](https://www.cs.ubc.ca/~andrones/param-training/). 

Here's a zipped version of my param file for convenience: 

[HKWS_rnasoft_param_repo.zip](doc_ASSETS/HKWS_rnasoft_param_repo.zip)

## LinearFold and LinearPartition  (`linear=True`)

Several arnie tools have the option to use LinearFold or LinearPartition to speed up calculation for long RNAs.  To set this up, download the following repos:

In `pfunc`, `mfe`, `bpps`, use the arg `linear=True`. This works for `package='vienna'` or `package='contrafold'`.

```
git clone https://github.com/LinearFold/LinearFold.git
cd LinearFold
make
cd ..
git clone https://github.com/LinearFold/LinearPartition.git
cd LinearPartition
make
```

This was most recently validated with LinearPartition commit `694b70a` and LinearFold commit `260c6bb`. Note that future commits that add more command-line options are likely to break this current wrapper. 

Check out these particular commits using

```
git clone https://github.com/LinearFold/LinearFold.git
cd LinearFold
git reset --hard 260c6bbb9bf8cc84b807fa7633b9cb731e639884

git clone https://github.com/LinearFold/LinearPartition.git
cd LinearPartition
git reset --hard 694b70a4d842af543d3f135af212fcb14a668bfd
```

add to arnie file:

```
# LinearFold build
linearfold: /path/to/LinearFold/bin

#LinearPartition build
linearpartition: /path/to/LinearPartition/bin
```

Note: EternaFold also includes a patch to use EternaFold parameters with LinearFold and LinearPartition. More information available with [EternaFold](https://eternagame.org/about/software).

## the TMP dir

Last but not least, you'll need to create a directory somewhere where your system calls to all of these packages will be deleting temporary files and then deleting them. Add to arnie file:

```
#directory to write temp files
TMP: /tmp
```
