# -​*- coding: utf-8 -*​-
# script to plot the results from R0s.dat
# no argument, everything is hardcoded
# @author: Laurent Hébert-Dufresne <lhebertd@uvm.edu>

# Packages
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import colorConverter as cc
import numpy as np

# =============================================================================
# Global parameters for the figures.
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 16
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Fira Sans", "PT Sans", "Open Sans", "Roboto", "DejaVu Sans", "Liberation Sans", "sans-serif"]
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.major.size"] = 8

#Load simulation results.
y = np.loadtxt('../data/SH_R0s.dat')

# colors
colors = ['#2b83ba', '#ba622b', '#80a339']

# Creates the figure and axis objects.
fig, ax = plt.subplots(1,1,figsize=(6,5.5), sharey=False)
fig.subplots_adjust(bottom=0.15)
#ax.set_xscale('log')
ax.set_xlim(0,6)
#ax.set_yscale('log')
ax.set_ylim(0,6)

#Show R0=1
x_coordinates = [0, 6]
y_coordinates = [1, 1]

ax.plot(x_coordinates, y_coordinates, c='Black', ls=':',  lw=1)

# Compute and plot final size
#    data = [data; baseR0, R0_NoSchools, R0_Lockdown, R0_LockdownWithXSchools, R0_LockdownWithSchools];
ax.plot(y[:,0], y[:,0], c='Black', ls='--',  lw=3, label=r'Baseline')
ax.plot(y[:,0], y[:,1], c=colors[1], ls='-',  lw=3, label=r'School closure alone')
ax.plot(y[:,0], y[:,2], c=colors[0], ls='-',  lw=3, label=r'Full lockdown')
ax.plot(y[:,0], y[:,5], c=colors[2], ls='-',  lw=3, label=r'Full school reopening')
ax.plot(y[:,0], y[:,3], c=colors[2], ls='-.',  lw=3,  label=r'Reopen $<$10y.o., partial above')
ax.plot(y[:,0], y[:,4], c=colors[2], ls='--',  lw=3,  label=r'Reopen $<$10y.o. only')

# Labels
ax.set_xlabel(r'Baseline $R_0$')
ax.set_ylabel(r'Effective $R_0$')

# Legend.
lines, labels = ax.get_legend_handles_labels()
ax.legend(lines, labels, loc='upper left', shadow=False, fancybox=False, prop={'size':13}, frameon=False, handlelength=3, numpoints=1)

# Save to file.
#plt.tight_layout(0.1)
fig.savefig("SH_Figure_R0s.png")
fig.savefig("SH_Figure_R0s.pdf")
