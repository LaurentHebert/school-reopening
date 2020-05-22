# LHD

Implementation of the model from Zhang et al. Science 2020 (DOI:10.1126/science.abb8001) to test different school reopening strategies.

### Codes

Integration of the ODEs is done in C++ using the Gnu Scientific Library, calculation of R0s in octave and plotting in python. 

The code tevol_source.cpp takes in a contact matrix and transmission rate and spits out a final epidemic size.

The octave code to test different strategies can probably be used in Matlab. There are no arguments, everything is set up as in the paper.

The python code to plot the figure plots the output of the octave code.

### Data

The relevant contact matrices are stored in data, are as temporary files. Population age-structure and susceptiblity are hardcoded in the codes.
