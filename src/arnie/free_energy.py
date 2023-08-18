import os, re, sys
import subprocess as sp
import random, string
import numpy as np
from .utils import *
from .pfunc import pfunc

DEBUG=False

# load package locations from yaml file, watch! global dict
package_locs = load_package_locations()

def free_energy(seq, constraint=None, package='vienna_2', T=37, coaxial=True, dna=False, beam_size=100,
		 pseudo=False, dangles=True, reweight=None, ensemble=True, param_file=None, linear=False,DEBUG=False):
	''' Compute free energy of RNA sequence. If structure is given, computes free energy of that structure. 
			Otherwise, returns MFE structure of sequence [NOT IMPLEMENTED YET].

		Args:
		seq (str): nucleic acid sequence
		constraint (str, optional): possible structure to constrain to in dot bracket notation
		T (float): temperature (Celsius), default 37

		ensemble (bool): to compute ensemble of constraint string or not.
			Just converts '.' to 'x' in string.
			If you want the free energy of just one structure,
			better practice is to use 'x' to denote unpaired. 


		motif (str): argument to vienna motif 
		beam_size (int): beam size for use in LinearPartition (Vienna, CONTRAfold, EternaFold only)
		dangles (bool): dangles or not, specifiable for vienna, nupack
                dna (bool): use SantaLucia model for DNA (NUPACK only)
		coaxial (bool): coaxial stacking or not, specifiable for rnastructure, vfold
		noncanonical(bool): include noncanonical pairs or not (for contrafold, RNAstructure (Cyclefold))
		pseudo (bool): include pseudoknot (nupack only)
		Implemented packages: 
		'vienna_1', 'vienna_2', 'contrafold'

		NB: doesn't multiply by kT for contrafold...
		
	Returns
		free energy (float)
	'''
	if not ensemble:
		constraint = constraint.replace('.','x')

	return pfunc(seq, package=package, T=T, dangles=dangles, coaxial=coaxial, pseudo=pseudo, dna=dna, beam_size = beam_size,
	 constraint=constraint, reweight=reweight, param_file=param_file, return_free_energy=True, linear=linear, DEBUG=DEBUG)

	# if package.lower().startswith('contrafold'):
	# 	Z_constrained = pfunc(seq, package=package, T=T, dangles=dangles, constraint=constraint,param_file=param_file)

	# 	return -1* np.log(Z_constrained) # .00198 is k in kcal/mol #0.0019899*(273+T) * 
	# else:
	# 	raise RuntimeError("%s `free_energy` not implemented yet" % package)
