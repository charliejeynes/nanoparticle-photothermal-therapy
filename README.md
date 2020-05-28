# nanoparticle-photothermal-therapy

This repository holds all code related to nanoparticle-photothermal-therapy currently underway at the University of Exeter. 

Monte Carlo simulations of light propagating through tissue is carried out in 'arc' (see pinned project from Freddy Wordingham on my git homepage)

Heat diffusion simulations, using input from 'arc', are found in the 'heat_diffusion.m' file

The use of the above is described in Jeynes et al. Biomolecules 2019, 9(8), 343; https://doi.org/10.3390/biom9080343

Nascent code for simulating diffusion of molecules in a microfludic geometry is found in the file 'diffusion_microfludic.py'
This is related to photothermal therapy in that the density of the nanoparticles found in a tumour is a dominating factor.
By simulating NP diffusion it may be possible to calculate the density of NPs in tumours and hence the photothermal effect. 
