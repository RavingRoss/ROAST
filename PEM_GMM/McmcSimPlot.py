# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 15:02:25 2022

@author: Steve
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import os
from os.path import exists
import statistics as st
import sys
import warnings
import astropy.stats as ast
from scipy import stats as sct
import random as rand
import time
import itertools as its

# Read the parameter file
# Home office desktop
param_dir = 'D:/Dropbox/ERAU/SARA Collaboration/Field roAp Photometry/Parameter_Files'
file_exists = exists(param_dir)
if file_exists == False:
    # School Office
    param_dir = 'E:/Dropbox/ERAU/SARA Collaboration\Field roAp Photometry/Parameter_Files'
    file_exists = exists(param_dir)
if file_exists == False:
    # Home Laptop
    param_dir = 'C:/Users/Steve/Dropbox/ERAU/SARA Collaboration/Field roAp Photometry/Parameter_Files'
    file_exists = exists(param_dir)
if file_exists == False:
    print ("Paramter directory not found.")
    sys.exit()

current_directory = os.getcwd()
os.chdir(param_dir)
param_file_list = os.listdir()
num_files = len(param_file_list)

print('  ')
print('Available Parameter Files')
print('--------------------------------------')
for i in range(0, num_files, 1):
    print(i, "  ", param_file_list[i])
print('--------------------------------------')
print('  ')
file_choice = input("Choose a parameter file [e.g. 0, 1, 2, 3..]: ")
file_choice = int(file_choice)
print('Using ' + param_file_list[file_choice])
param_file = param_file_list[file_choice]
param_path = param_dir + "/" + param_file

with open(param_path) as par:
    for line in par:
        if line.strip():
            word = line.split(':\t')
            word1 = word[0]
            if word1[0] != '#':
                word2 = ((word[1].split('\n'))[0]).split('\t')
                n2 = len(word2)
                word2 = word2[n2-1]

                # Analysis case keywords
                if word1 == 'TESS_path':
                    path = word2
                if word1 == 'Data_path_1':
                    path = word2
                    file_exists = exists(path)
                if word1 == 'Data_path_2' and file_exists == False:
                    path = word2  
                    file_exists = exists(path)
                if word1 == 'Data_path_3' and file_exists == False:
                    path = word2            
                    file_exists = exists(path)  
                if word1 == 'target':
                    target = word2
                if word1 == 'sector':
                    sector = word2
                if word1 == 'case':
                    case = word2
                if word1 == 'subcase':
                    subcase = word2
                if word1 == 'itermax':
                    itermax = int(word2)
                if word1 == 'kauth':
                    kauth = int(word2)
                if word1 == 'tol':
                    tol = float(word2)    
                if word1 == 'h':
                    h = float(word2)    
                if word1 == 'a':
                    a = float(word2)  

if file_exists == False:
    print("McMc simulation directory not found.")
    sys.exit()
    
os.chdir(current_directory)

# Analysis case path
if target == '':
    target = input("Target [e.g. KIC7582608]: ")
if sector == '':
    sector = input("Sector [e.g. SEctor 14 2019]: ")
if case == '':
    case = input("Analysis case [e.g. Case1, Case 2, ...]: ")
if subcase == '':
    subcase = input("Subcase [e.g. /Days/, /Months/]: ")

file = input(
    "Montecarlo sim results file [e.g.KIC7582608_BJD2459700_Mcmc_log.txt ]: ")

directory = "".join([path, target, "/", sector, "/", case, subcase])
file_path = "".join([directory, file])

current_directory = os.getcwd()
os.chdir(directory)

# Count lines in Montecarlo sim file
# There are *hlines* header lines plus *nentries* lines of data
with open(file_path, 'r') as fp:
    hlines = 0
    while (line[0] != " Proc"):
        line = fp.readline().split("#")
        hlines += 1
        if (line[0].split(":"))[0] == "      (successfully completed)":
            nentries = int(((line[0].split(":"))[1].split("\n"))[0])

data = np.loadtxt(file_path,  dtype='float', comments='#', delimiter=None,
                  converters=None, skiprows=hlines, usecols=None, unpack=False, ndmin=0, max_rows=nentries)
rows, cols = data.shape
nfreq = int((cols - 2)/3)
data = np.transpose(data)

# Global varaibles
ndata = rows
nbins = int(200*ndata/10000) # Number of histogram bins
print ("Ndata = ", ndata)
print("Nbins = ", nbins)

def looksee(x, quantity, units):
    # Plot the data PDF to identify features 
    lower = min(x)
    upper = max(x)
    
    prange = np.linspace(lower, upper, nbins)
    count, bins, bars = plt.hist(x, nbins, density=True, log=True, color='blue', alpha=0.7)
    x_label = units
    y_label = "Log Probability Density(in Frequency" + units + ")"
    ptitle = "Look-See Plot of  " + quantity    
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.suptitle(ptitle, y=0.92, fontsize=14)
    plt.show()
    
    done = "No"
    while (done == "No"):
        fig = plt.figure(figsize = (11,8))
        ax0 = fig.add_subplot(111)
        count, bins, bars = plt.hist(x, nbins, density=True, log=True, color='blue', alpha=0.25)
        x_label = units
        y_label = "Log Probability Density(in Frequency/mhz)"
        ptitle = "Look-See Plot of  " + quantity    
        plt.xlabel(x_label, fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        plt.suptitle(ptitle, y=0.92, fontsize=14)
    
        parameters = []
        proportions = []
        parameters = [float(item) for item in input("Trial Mean(s) & Stdev(s) of " + quantity + " [e.g. 2.5 0.005 2.15 0.005]: ").split()]
    
        ngauss = int(len(parameters)/2)
    
        j=0
        for i, c in zip(range(ngauss), ['r', 'g', 'b', 'c', 'm', 'y']):
            xl = parameters[j] - 5*parameters[j+1]
            xu = parameters[j] + 5*parameters[j+1]
            plt.axvline(x=xl, ymin=1e-16, linewidth=2, color=c)
            plt.axvline(x=xu, ymin=1e-16, linewidth=2, color=c) 
            j = j + 2
        plt.show()
        done = input("Are you satisified with the placements? [Yes or No] : ")  
    
    proportion = 1/ngauss
    proportions = [0]*ngauss
    for i in range(0, ngauss):
        proportions[i] = proportion
    print ("Initial Proportions =", proportions)
    
    return parameters, proportions, nbins, lower, upper, ndata

def logLikelyhood(x, ngauss, pi, gdev):
    # Calculate the log(likelyhood)
    likelyhoods = [0]*ngauss
    for c in range(ngauss):
        likelyhoods[c] = pi[c]*(gdev[c].pdf(x))+1e-10
    likelyhoods = np.nan_to_num(likelyhoods, nan=1e-10)
    log_likelyhood = np.sum(np.log(np.sum(np.array(likelyhoods), axis=0)))
    return log_likelyhood

def confidence_interval(N, S, level, DoF):
    # N = number of data used to calc. sigma
    # S is the sample mean
    # level is the confidence level, 0 <= level <= 1
    # DoF = degrees of freedom
    t = sct.t.interval(level, DoF)[1]
    confint = (t*S/np.sqrt(N))
    
    return confint

def bin_errors(x, nsamples, sample_size):
    # Estimating the uncertainties on the bin values
    # nsamples samples of size sample_size are taken from the list of 
    # random deviates x.
    # A histogram is constructed of each sample and the standard deviations of
    # of each bin calcuated.
    nxdata = int(len(x)/nsamples)
    nruns = int(sample_size/nsamples)
    
    esigmas = []
    for j in range(0, nruns):
        # Resampling the data. nsample random samples.
        xsample = []
        pdxsample = []
        
        for i in range (0, nsamples):
            xsample = np.append(xsample, rand.sample(list(x), nxdata))
            pd, bin_edges = np.histogram(xsample, nbins, density=True)
            pdxsample = np.append(pdxsample, pd)
        
        pdxsample = np.array(np.reshape(pdxsample, [nsamples, nbins]))
        esigmas = np.append(esigmas, np.std(pdxsample, axis=0, dtype=None, \
                            out=None, ddof=1))

    esigmas = np.reshape(esigmas, [nruns, nbins])    
    esigmas = np.mean(esigmas, axis=0)
    esigmas = confidence_interval(nsamples, esigmas, 0.68, nsamples-2)
    return esigmas

def EM(x, proportions, gaussians, quantity, units, kauth, esigmas):
    #
    #=============================================================
    # Expectation - (Log Likelyhood) Maximization (EM) algorithm
    #=============================================================
    #
    # Algorith taken from
    # "Expectation maximization and Gaussian mixture models",
    # By Tobias Schlagenhauf. Last modified: 22 Feb 2022.
    # URL: https://python-course.eu/machine-learning/expectation-maximization-and-gaussian-mixture-models-gmm.php
    #
    #--------------------
    # Expectation step
    #--------------------
    #
    # x is the data array 
    # gaussians is a list of mean and sigma pairs 
    # specifying the starting guesses for the Gaussians
    # e.g. [mu0, sigma0, mu1, sigma1, mu2, sigma2]
    # proportions is an array holding the fraction of the distribution
    # belonging to eaxch gaussian.
    # kauth is the maximum number of iterations allowed.
    # Clip the data per the upper and lower limits below
    
    pis = [0]*kauth
    mus = [0]*kauth
    sigs = [0]*kauth
    logLs = [0]*kauth
    
    lower = min(x)
    upper = max(x)
    
    prange = np.linspace(lower, upper, nbins)
    xlower = lower
    xupper = upper
    
    ngauss = int(len(gaussians)/2)
    r = np.zeros((ndata, ngauss))
    epsilon = 1000
    log_likelyhood = -1e12
    
    # Initialize the ngauss gaussian deviates.
    j=0
    gdev = [0]*ngauss
    mu_c = [0]*ngauss
    sigma_c = [0]*ngauss
    for i in range(0, ngauss):
        mu_c[i] = gaussians[j]
        sigma_c[i] = gaussians[j+1]
        gdev[i] = sct.norm(loc=mu_c[i], scale=sigma_c[i])
        j = j + 2
    pi = np.array(proportions)
    must = mu_c
    sigst = sigma_c
    index_max_pi = np.argmax(np.array(pi))
    
    
    # Iterate
    k = 1 

    while (epsilon > tol and k <= kauth):
        old_log_likelyhood = log_likelyhood
        
        # Limit the PDFs to this range
        prange = np.linspace(lower, upper, nbins)
        
        # Do not let the proportions go negative
        if any(x < 0 for x in pi):
            index_max_pi = np.argmax(np.array(pi))
            pi = [0]*ngauss
            pi[index_max_pi] = 1.
            
        #Probability that each datapoint belongs to each gaussian distribution
        for c,g,p in zip(range(ngauss), gdev, pi):
            r[:, c] = p*g.pdf(x)
    
        # Normalize the probabilities so that each row sums to one. Then weight it
        # by the fraction of points belonging to cluster (gaussian) c.
        for i in range(len(r)):
            rowsum = np.nansum(r, axis=1)[i]
            if rowsum > 0:
                r[i] = r[i]/(np.nansum(pi)*rowsum)
        fig = plt.figure(figsize = (11,8))
        ax0 = fig.add_subplot(111)
        
        Gpdf =[0]
        gpdf = [0]*ngauss
        for i in range(ngauss):
            gpdf[i] = gdev[i].pdf(prange)
            Gpdf = Gpdf + gpdf[i]*pi[i]  # Sum the histograms
         
        # Plot the data PDF    
        binvalues, binedges, patches = plt.hist(x, nbins, density=True, color='blue', alpha=0.25)
        
        # Plot the mixture model PDFs
        for g, p, c in zip(gpdf, pi, ['r', 'g', 'b', 'c', 'm', 'y']):
            ax0.plot(prange,g*p,c=c,zorder=0)
        
        ax0.plot(prange,Gpdf,color="k",zorder=0)
        
        binwidth = binedges[1] - binedges[0]
        scaleFactor = binwidth*ndata
        
        # Calculate the goodness of fit
        Nobs = binvalues*scaleFactor # number per bin
    
        Nmodel = Gpdf*scaleFactor
        #Nmodel = gpdf[index_max_pi]*scaleFactor*pi[index_max_pi]
        Tmodel = sum(Nmodel)
        
        # Estimating the uncertainties on the bin values of the model
        bin_variances = np.square(pi[0]*bin_errors(gdev[0].rvs(ndata), ns, ndata/ns))
        for i in range(1, ngauss):
            bin_variances = np.add( bin_variances, np.square(pi[i]*bin_errors(gdev[i].rvs(ndata), ns, ndata/ns) ) )
        model_errors = np.sqrt(bin_variances)
        
        sigmas = np.sqrt(esigmas**2 + model_errors**2)*scaleFactor # estimated in the initialization step.
        GoF = chiSquaredPerDOF(Nobs, Nmodel, sigmas, nbins-2*ngauss)
        
        # Histogram Labels
        binwidth_string = "Bin width = " + \
            str(format(binwidth, ".4E"), ) + " " + units
        #compno_string = "Component "+ str(index_max_pi) + "\n"
        compno_string = "Sum of " + str(ngauss) + " Components\n"
        title_string = compno_string + "N = " + str(int(Tmodel)) + "/" + str(ndata) + "\n"
        title_string = title_string + r'$\chi^2_\nu$' + "= " + str(round(GoF, 3)) \
            + "\n" + binwidth_string  
        
        legends = [0]*(ngauss+2)
        for i in range(0, ngauss):
            legends[i] = "Component:" + str(i) + "\n" + r'$\mu$' + " = " + str(round(mu_c[i], 5)) + \
            " " + units + "\n" + r'$\sigma$' + " = " + \
            str(round(sigma_c[i], 5)) + " " + units + "\n" \
            + r'$\pi$' + " = " + str(round(pi[i], 5))
        legends[ngauss] = "Sum of components"     
        legends[ngauss+1] = "Mcmc-simulated" + " " + quantity
        plt.title(title_string, x=0.01, y=0.83, fontsize=12, loc="left")
        plt.legend(legends, fontsize=12, loc="upper right")
        
        x_label = quantity + " [" + units + "]"
        y_label = "Probability Density" + "(Frequency/" + units + ")"
        ptitle = "Comparison of Simulated " + quantity + " and a Gaussian Mixture Model "   
        plt.xlabel(x_label, fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        plt.suptitle(ptitle, y=0.92, fontsize=14)
        plt.axis([xlower, xupper, None, None])  
        
        fignm = quantity + "Fit.jpeg"
        filename = directory + "/" + fignm
        # fig = plt.gcf()
        plt.savefig(filename, dpi=600, bbox_inches='tight')
        
        plt.show()
        
        #------------------------
        # Maximization step
        #-----------------------
        #
        #
        if any(x == 0 for x in pi):
            # Save the last kauth iterations
            pis[k-1] = np.array(pi)
            mus[k-1] = mu_c
            sigs[k-1] = sigma_c
            logLs[k-1] = log_likelyhood
            break
        
        # The number of points belonging to each gaussian distribution.
        m_c = []
        for c in range(len(r[0])):
            m = np.nansum(r[:, c])
            if m == 0:
                m = 1e-16
            m_c.append(m)
        
        # The fraction of points belonging to each gaussian distribution.    
        pi_c = []
        for m in m_c:
            pi_c.append(m/np.nansum(m_c))
    
        # Convert the row vector into a column vector
        x_T = x.reshape(ndata, 1) 
    
        # Update the means
        mu_c = np.nansum(x_T*r, axis=0)/m_c
        
        # Calculate deviations
        dev = x_T-mu_c[c]
        
        # Update the sigmas
        var_c = []
        for c in range(len(r[0])):
            r_c_T = np.array(r[:, c]).reshape(ndata,1)
            var_c.append((1/m_c[c])*np.dot((r_c_T*dev).T, dev) )
            
        sigma_c = np.sqrt(var_c).flatten()
         
        # Update the gaussian distributions
        gdev = [0]*ngauss
        for i in range(0, ngauss):
            mu = mu_c[i]
            sigma = sigma_c[i]
            gdev[i] = sct.norm(loc=mu, scale=sigma)
    
        # Update the proportions
        pi = pi_c
        
        # Calculate the log likileyhood
        log_likelyhood = logLikelyhood(x, ngauss, pi, gdev)
        
        print("Iteration: ", k, "/", kauth, " Mu: ", mu_c, " Sig: ", sigma_c, \
          "log(L) = ", log_likelyhood)
        
        # Save the last kauth iterations
        pis[k-1] = np.array(pi_c)
        mus[k-1] = mu_c
        sigs[k-1] = sigma_c
        logLs[k-1] = log_likelyhood
        
        epsilon = abs(log_likelyhood - old_log_likelyhood )
        
        # Adjust the x-axis limits
        xlower = mu_c[index_max_pi] - 7*sigma_c[index_max_pi]
        xupper = mu_c[index_max_pi] + 7*sigma_c[index_max_pi]
        k = k + 1
      
    # Niters-dimensional point in solution space
    niters = k - 1
    if niters == 0:
        pist = [0]*(1)
        must = [0]*(1)
        sigst = [0]*(1)
        #logLst = [0]*(1)
        pist[0] = np.array(pi)
        must[0] = np.array(mu_c)
        sigst[0] = np.array(sigma_c)
        
        logLst = logLs
        axis = 1
    else:
        pist = [0]*(niters)
        must = [0]*(niters)
        sigst = [0]*(niters)
        logLst = [0]*(niters)
        for i in range(0,niters):
            pist[i] = pis[i]
            must[i] = mus[i]
            sigst[i] = sigs[i]
            logLst[i] = logLs[i]
        axis = 1
    
    P = np.concatenate((pist, must, sigst), axis=axis)
    
    return P, logLs, niters, GoF

def ChauvenetCrit(x):
    outliers = []
    clean = []
    indices = []
    N = len(x)
    mu = np.nanmean(x)
    sigma = np.nanstd(x)
    for i in range(0, N):
        Z = abs(x[i]-mu)/sigma
        test = N*math.erfc(Z/np.sqrt(2))
        if test < 0.5:
            outliers = np.append(outliers, x[i])
        else:
            clean = np.append(clean, x[i])
            indices = np.append(indices, i)
    return clean, outliers, indices


def chiSquaredPerDOF(o, c, sigmas, DOF):
    ocln = []
    ccln = []
    sigcln = []
    index = np.nonzero(sigmas)
    for i in index:
        ocln = np.append(ocln, o[i])
        ccln = np.append(ccln, c[i])
        sigcln = np.append(sigcln, sigmas[i])

    r = np.divide(ocln - ccln, sigcln)
    chiSquared = np.sum(r**2)
    CSPD = np.divide(chiSquared, DOF)
    
    return CSPD
  
def Plot3d(x, y, z, ptitle, x_label, y_label, z_label, figno):

    figtitle = 'fig' + str(figno)
    figtitle_out = target + "_Mcmc_Plot_" + figtitle + '.png'
    figtitle = plt.figure() 
        
    ax = figtitle.add_subplot(111, projection = "3d")
    ax.scatter(x, y, z, s=1)

    legend3 = "Mcmc Estimates"

    if ptitle != "":
        ax.set_title(ptitle)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    plt.legend([legend3], fontsize=10, loc="upper right")
    if figno != "":
        figtitle.savefig(figtitle_out, dpi=600, bbox_inches='tight')
        
        
def SaveResults(McmcStats, directory, quantity, units, freqno, ngauss):
        header1 = "F   Proportion      " + quantity + "           " +  quantity + ".Sigma      GoF"
        header2 = "                   [" + units + "]            [" + units + "]"
        results_out = directory + 'Mcmc_Statistics.txt'
        print ("Saving resuls to ", results_out)
        f = open(results_out, 'a')
        f.write("MCMC STATISTICS\n")
        f.write("  \n")
        f.write(header1 + "\n")
        f.write(header2 + "\n")

        print(" ")
        print("MCMC STATISICS ")
        
        print(header1)
        print(header2)
        for k in range(0, ngauss*4, 4):
            print("%2d   %4.5f  %3.14f  %3.14f  %4.3f" % (freqno, McmcStats[k], McmcStats[k+1], McmcStats[k+2], McmcStats[k+3]))
            f.write("%2d   %4.5f  %3.14f  %3.14f  %4.3f \n" % (freqno, McmcStats[k], McmcStats[k+1], McmcStats[k+2], McmcStats[k+3]))
        f.close()

def Parabolic_Accelerator(data, quantity, units):
        #---------------------
        # PEM starting run
        #---------------------
        #
        # Taken from 
        # "Parabolic acceleration of the EM algorithm"
        # A. Berlinet and C. Roland, Published Online 16 May 2008 
        # in Stat Comput (2009) 19: 35–47 
        # DOI 10.1007/s11222-008-9067-x
        #
        McmcStats = []
        iterations = 0
        ngauss = int(len(gaussians)/2)
        
        # Kauth is set in the paramter file.
        P, L, niters, GoF = EM(data, proportions, gaussians, quantity, units, kauth, esigmas) 
        iterations = iterations + niters
        # Control Points and log(L) of most recent but one update
        CP = [0]*3
        for i in range(0,3):
           j = niters - 3 + i
           CP[i] = P[j] # Control points
            
        P0 = np.array(CP[0])
        P1 = np.array(CP[1])
        P2 = np.array(CP[2])
        
        # Already converged so don't run the PEM.
        nps = len(P)-1
        if np.linalg.norm(P2-P1)**2 <= tol**2:
            for i in range(0, ngauss):
                McmcStats = np.append(McmcStats, [P[nps][i], P[nps][ngauss+i], P[nps][2*ngauss+i], GoF])
            print ("Parabolic acceleration not applied as the solution has already converged.")
            return McmcStats
            
        L2_old = -1e200    
        
        ngauss = int(len(P2)/3) # There are three parameters to optimize.
        
        gdev = [0]*ngauss
        pi = [0]*ngauss
        mu = [0]*ngauss
        sigma = [0]*ngauss
        mu_sigma = [0]*ngauss
        
        for iter in range(1, itermax):
            
            for i in range(0, ngauss):
                pi[i] = P2[i]
                mu[i] = P2[i+ngauss]
                sigma[i] = P2[i+2*ngauss]
                gdev[i] = sct.norm(loc=mu[i], scale=sigma[i])    
        
            L2 = logLikelyhood(data, ngauss, pi, gdev)
            
            # Geometric search
            # Tunable search parameters
            # h = Search grid mesh size
            # a = Ratio
            
            # Search step
            i = 1
            t = 1 + h * a**i 
        
            # Exploration step
            Pnew = P2*t**2 + 2*t*(1-t)*P1 +  P0*(1-t)**2
            
            # Unpack the Point Pnew to get the updated pis, mus, and sigmas
            # So as to make the updated gaussian deviates
             
            for i in range(0, ngauss):
                pi[i] = Pnew[i]
                mu[i] = Pnew[i+ngauss]
                sigma[i] = Pnew[i+2*ngauss]
                gdev[i] = sct.norm(loc=mu[i], scale=sigma[i])    
            
            Lnew = logLikelyhood(data, ngauss, pi, gdev)
            
            if Lnew <= L2:
                print("Lnew <= L2")
                P0 = P2
            
                for i in range(0, ngauss):
                    pi[i] = P0[i]
                    mu[i] = P0[i+ngauss]
                    sigma[i] = P0[i+2*ngauss]
                    mu_sigma[i] = [mu[i], sigma[i]]
                
                musigma = list(its.chain.from_iterable(mu_sigma)) # Flatten the parameter list
                P1, L, niters, GoF = EM(data, pi, musigma, quantity, units, 1, esigmas) 
                iterations = iterations + niters
                P1 = np.array(list(its.chain.from_iterable(P1))) # Flatten the solution point
            
                for i in range(0, ngauss):
                    pi[i] = P1[i]
                    mu[i] = P1[i+ngauss]
                    sigma[i] = P1[i+2*ngauss]
                    mu_sigma[i] = [mu[i], sigma[i]]  
                
                musigma = list(its.chain.from_iterable(mu_sigma)) # Flatten the parameter lists    
                P2, L, niters, GoF = EM(data,  pi, musigma, quantity, units, 1, esigmas) 
                Pout = P2
                nps = len(Pout)-1
                iterations = iterations + niters
                P2 = np.array(list(its.chain.from_iterable(P2))) # Flatten the solution point
            else:
                
                while Lnew > L2:
                    print ("Lnew > L2")
                    Pold = Pnew
                    L2 = Lnew
                    i = i + 1
                    t = 1 + h * a**i 
                    Pnew = P2*t**2 + 2*t*(1-t)*P1 +  P0*(1-t)**2
                    
                    for i in range(0, ngauss):
                        pi[i] = Pnew[i]
                        mu[i] = Pnew[i+ngauss]
                        sigma[i] = Pnew[i+2*ngauss]
                        mu_sigma[i] = [mu[i], sigma[i]]  
                        gdev[i] = sct.norm(loc=mu[i], scale=sigma[i])    
                
                    Lnew = logLikelyhood(data, ngauss, pi, gdev)    
                
                P0 = P1
                P1 = P2
                
                musigma = list(its.chain.from_iterable(mu_sigma)) # Flatten the parameter lists  
                P2, L, niters, GoF = EM(data,  pi, musigma, quantity, units, 2, esigmas) 
                
                nps = len(P2)-1
                if nps >= 0:
                    Pout = P2
                    P2 = P2[nps]
                else:
                    Pout = P1
                P2 = np.array(P2)   
                iterations = iterations + niters
            
            print("Lnew = ", Lnew)
            if (abs(L2-L2_old) <= tol or nps == -1):
                print ("Total Iterations =",  iterations)
                for i in range(0, ngauss):
                    McmcStats = np.append(McmcStats, [Pout[nps][i], Pout[nps][ngauss+i], Pout[nps][2*ngauss+i], GoF])
                
                return McmcStats
            else:
                L2_old = L2
        
        print ("Total Iterations =",  iterations)
        for i in range(0, ngauss):
            McmcStats = np.append(McmcStats, [Pout[nps][i], Pout[nps][ngauss+i], Pout[nps][2*ngauss+i], GoF])   
        return McmcStats
    
ns = 10 # number of subsamples
    
#amplitudes, outliers, indices = ChauvenetCrit(amplitudes)

labels = [""]*3
labels[0] = ["Frequencies", "mHz"] 
labels[1] = ["Amplitudes", "mmag"] 
labels[2] = ["Phases", "x2" + r'$\pi$' + " radians"] 

for f in range(0, nfreq*3, 3):
    freqno = int(1 + f/3)
    amplitudes = np.multiply(data[f+3], 1000)
    Plot3d(data[f+2], data[f+4], amplitudes, "F" + str(freqno) + " Amplitude-Frequency-Phase Space", \
         labels[0][0] + "[" + labels[0][1] + "]", \
         labels[2][0] + "[" + labels[2][1] + "]", \
         labels[1][0] + "[" + labels[1][1] + "]", "")
    plt.show()
    
    for d,l in zip([data[f+2], amplitudes, data[f+4]], labels):
        # Look at the parameter space to estimate number of gaussians needed
        gaussians, proportions, nbins, lower, upper, ndata = looksee(d, l[0], l[1])
    
        # Estimating the uncertainties on the bin values 
        esigmas = bin_errors(d, ns, ndata/ns)
    
        # PEM Algorithm
        McmcStats = Parabolic_Accelerator(d, l[0], l[1])
        
        # Save means and sigmas to a file
        ngauss = int(len(gaussians)/2)
        SaveResults(McmcStats, directory, l[0], l[1], freqno, ngauss)
        
sys.exit()



