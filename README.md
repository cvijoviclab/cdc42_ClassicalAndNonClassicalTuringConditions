# About 
This is the Github repositry of the manuscript ``Cell polarisation in a bulk-surface model can be driven by both classic and non-classic Turing instability''[1]. The model simulates the process called Cdc42-mediated cell polarisation which is essential for multiple processes such as cell division and cell migration for a wide variety of cells ranging from yeast to human beings. During this process, the protein Cdc42 belonging to the class of Rho-GTPases gets activated and inactivated as well as being transported from the interior of the cell (the cytosol) to the surface of the cell (the cell membrane) through a set of chemical reactions. Simultaneously,it moves in space through diffusion, and it is this combination of reaction and diffusion that ultimately causes the aggregation of the active component of cdc42 in one spot on the cell membrane called the pole. Recent attempt to model this process has included a three dimensional description of the cell including both the cytosol (i.e. a bulk) and the membrane (i.e. a surface) [2,3], and for this type of **bulk-surface models** there are two types of mechanisms that can cause pattern formation. These are so called *classic Turing instability* (Fig 1A) and a more recent version referred to as *non-classic Turing instability* (Fig 1B).  





![**Figure 1:** The evolution of a pole being a single spot of active cdc42 for two different cases: (A) Classic Turing instability and (B) Non-classic Turing instability.](./Figures/evolutionPattern/evolutionPattern.png "Title")


In this paper, based on the structure of a previous model [2,3] we propose a simpler and (what we argue for) a more biologically realistic bulk surface model of Cdc42-mediated cell polarisation. We conduct mathematical analysis in order to conclude that the solutions to the proposed model are physically reasonable. Then, we conduct a numerical investigation of the parameter space in order to map out the sets of parameters that give rise to the classic and non-classic type of Turing instability. Lastly, based on this numerical description of the parameter space, we select parameters giving rise to cell polarisation (at least theoretically) and we validate our results using numerical simulations of the reaction diffusion model. The numerical simulations are based on a combined approach using the *Finite Element Method (FEM)* [4] in space together with an *adaptive Finite Difference (FD)* [5] approach in time. Using these simulations, we investigate the effect of changing the involved rate parameters on three different properties: the minimal and maximal concentration of cdc42 on the cell membrane, the polarisation time being the time it takes for the pole to be formed and the ratio between the area of the pole relative to the surface of the cell membrane. In particular, we investigate two key parameters being the relative diffusion between the inavtive and active state of cdc42 denoted d and the relative reaction strength gamma which is directly proportional to the surface area of the cell membrane.  



In this repository, a script for generating the numerical description of the parameter space as well as two scripts for simulating the cdc42-mediated cell polarisation caused by diffusion driven instability are available in the ``Code'' folder. Here, code for processing and formatting the meshes, i.e. the spatial discretisations, used in the FEM are also found. In addition, code for analysing the results, i.e. the output data, from the FEM-FD simulations is presented as well as code for generating multiple figures presented in the article. 

## Repository structure
The following folders are contained in the repositry: 

- **Code:** The folder contains the scripts for running the FEM-FD algorithm for simulating the evolution of the RD <br>
system, it contains the script for generating a plot of the parameter space which determines which parameters that <br>
give rise to classic and non-classic Turing instability, it contains the script for creating the meshes used in the <br>
FEM-FD algorithm which corresponds to the spatial approximation of the cell, it contains a script for calculating the  <br>
steady states of the homogeneous system which are used when setting the initial conditions in the FEM-FD algorithm and <br>
lastly it contains a script for analysing the data obtained from a running the ``increasing gamma'' and ``increasing d''<br>
experimental designs. 
- **Results:**
- **Figures:**

Each of the sub folders in these three major folders contains their own ``README.md'' files which explains the details of how each script should be executed. 

## Requirements for reproducing the result



### Operating system


### Programming languages



### Software






















# References
1. Borgqvist, Johannes, et al. "Cell polarisation in a bulk-surface model can be driven by both classic and non-classic Turing instability." bioRxiv (2020).
2. Rätz, Andreas, and Matthias Röger. "Turing instabilities in a mathematical model for signaling networks." Journal of mathematical biology 65.6-7 (2012): 1215-1244.
3. Rätz, Andreas, and Matthias Röger. "Symmetry breaking in a bulk–surface reaction–diffusion model for signalling networks." Nonlinearity 27.8 (2014): 1805.
4. Larsson, Stig, and Vidar Thomée. Partial differential equations with numerical methods. Vol. 45. Springer Science & Business Media, 2008.
5. Cheney, E. Ward, and David R. Kincaid. Numerical mathematics and computing. Cengage Learning, 2012.
