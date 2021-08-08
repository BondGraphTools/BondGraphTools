BondGraphsTools from 10,000 feet
================================

BondGraphTools is designed to be a modelling tool, to allow 
the user to programmatically describe the power storing, transforming and dissipating elements of a system, and generate simplified/reduced symbolic equations for the system as a whole.


Model Capture
-------------
A bond graph model is created by:
 1. Instantiating a new bond graph - a composite container for the system, or part of a larger system, that the user wishes to model.
 2. Existing components (either atomics drawn from standard libraries, or composites built via this process) are added to the bond graph model, describing the various processes (eg, dissipation, or energy storage) active within the system in question.
 3. Connections are defined between the ports of components, capturing the channels through which power flows. 
 4. If necessary, ports are `exposed`, defining a power port outside the container which mirrors the power flowing through a corresponding interior port.

Once a model is created, it is ready for reduction and/or simulation, or for integration with other models to form a larger, more compex system. 


Model Reduction
---------------
In the context of BondGraphTools, model reduction refers to the process by which the constituative relations (equations which constrain power flows and storage within a given component) and connectivity constraints (ie, bonds) are passed through a series of symbolic algebraic simplifications so as to generate a lower-dimensional represenation of the system.

The model capture procedure results in composite state $X$ containing all dynamic variables (eg, position/momentum), control variables and power variables (efforts and flows). The behaviour of the system is then governed by the level sets of $\Phi(X) = 0$. By exploiting inherient structures (such a sparsity and a only small amounts of nonlinearly), an invertible co-ordinate transform $Y = P(X)$ can be found such that $dim(Y) \ll dim(X)$, along with a reduced order constitutive relation for the composed system $\tilde{\Phi}(Y) = 0$, where $\tilde{\Phi}$ is 'simpler' in the sense that it contains a minimal set of relations describing the systems dynamics.


This model reduction happens behind the scenes when the user asks for the constiutive relations of their newly-constructed model, or in prepation for running numerical simulation.
