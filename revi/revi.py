import numpy as np

from glob import glob
import os, sys, pickle, requests

from arnie.pfunc import pfunc
from arnie.bpps import bpps
from arnie.sample_structures import sample_structures
from arnie.utils import convert_dotbracket_to_matrix
from decimal import Decimal
from arnie.free_energy import free_energy
from copy import copy

class REVI:
    '''RNA Ensemble Variational Inference.
    Input: sequence (str): RNA sequence
        n_samples (int): number of samples stochastically sampled each iteration.
        n_iters (int): number of iterations
    Returns:
    '''
    
    def __init__(self, seq, n_samples=10, n_iters=10, kT=0.6):
        self.seq = seq
        self.n_samples = n_samples
        self.n_iters = n_iters
        self.kT = 0.6
        self.scale = 1
        self.alpha_vectors = []
        self.perturbation_vectors = [np.zeros(len(seq))]
        self.grad_omega = [np.zeros(len(seq))]
        self.init_dot_plot = []
        self.final_dot_plot = []
        self.dG_v_uncon = []
        self.dG_c_uncon = []
        self.dG_corr = []
        self.posteriors = []
        
    def alpha_mean_(self):
        return np.mean(self.alpha_vectors, axis=1)
        
    def write_perturbation_file_(self, iter_ind):
        with open("perturb.%d.dat" % iter_ind,'w') as f:
            for (i, dE) in enumerate(self.perturbation_vectors[-1]):
                f.write("E %d 0 1 %.9f\n" % (i+1, dE))

    def reweight(self):
        print("REVI is working:")
        for iter_ind in range(self.n_iters):
            #print("iter %d" % iter_ind)
            if iter_ind==0:
                perturb_file=None
            else:
                perturb_file = 'perturb.%d.dat' % (iter_ind-1)
        
            self.reweight_step_(reweighting_file = perturb_file)
            self.write_perturbation_file_(iter_ind)
                
    def reweight_step_(self,reweighting_file=None):

        sampled_structs_uncon = sample_structures(self.seq, n_samples=self.n_samples, package='vienna_2',
                                                  nonredundant=True, reweight=reweighting_file)
        
        dG_v_uncon = np.array([free_energy(self.seq, constraint = s, package='vienna_2',
                                           reweight=reweighting_file) for s in sampled_structs_uncon])
        
        dG_c_uncon = np.array([free_energy(self.seq, constraint = s, package='contrafold') for s in sampled_structs_uncon])

        Z_v_uncon = np.sum(np.exp(-1/self.kT*dG_v_uncon))
        Z_c_uncon = np.sum(np.exp(-1*dG_c_uncon))

        Z_v_con, Z_c_con = [],[]

        for i in range(len(self.seq)):
            constraint = ['.']*len(self.seq)
            constraint[i] = 'x'
            constraint = ''.join(constraint)

            sampled_structs = sample_structures(self.seq, n_samples=self.n_samples, package='vienna_2',
                                                constraint=constraint, nonredundant=True, reweight=reweighting_file)

            dG_v = np.array([free_energy(self.seq, constraint = s, package='vienna_2',
                                         reweight=reweighting_file) for s in sampled_structs])
            
            dG_c = np.array([free_energy(self.seq, constraint = s, package='contrafold') for s in sampled_structs])

            Z_v_con.append(np.sum(np.exp(-1/self.kT*dG_v)))
            Z_c_con.append(np.sum(np.exp(-1*dG_c)))
            
        err = Z_c_con/Z_c_uncon - Z_v_con/Z_v_uncon
        self.dG_v_uncon.append(dG_v_uncon)
        self.dG_c_uncon.append(dG_c_uncon)
        self.dG_corr.append(np.corrcoef(dG_v_uncon, dG_c_uncon)[0][1])

        self.grad_omega.append(err)
        
        perturbation_vector = self.perturbation_vectors[-1] - self.scale*self.grad_omega[-1]
    
        self.perturbation_vectors.append(perturbation_vector)
        self.posteriors.append(Z_v_con/Z_v_uncon)

    def plot_perturbation_vectors(self):
        if len(self.perturbation_vectors) == 0:
            raise RuntimeError("Not estimated! Can't plot perturbation vectors")
        else:
            colors=sns.color_palette('magma_r',estimator.n_iters)
            for i,vec in enumerate(estimator.perturbation_vectors[1:]):
                plot(vec,c=colors[i],label=i)
                
    def plot_grad_omega(self):
        if len(self.grad_omega) == 0:
            raise RuntimeError("Not estimated! Can't plot posteriors")
        else:
            colors=sns.color_palette('magma_r',self.n_iters)
            for i in range(self.n_iters):
                plot(self.grad_omega[i],c=colors[i],label=i)   
                
    def plot_posteriors(self):
        if len(self.posteriors) == 0:
            raise RuntimeError("Not estimated! Can't plot posteriors")
        else:
            colors=sns.color_palette('magma_r',self.n_iters)
            for i,vec in enumerate(self.posteriors):
                plot(vec,c=colors[i],label=i)