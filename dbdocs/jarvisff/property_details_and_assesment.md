Using jarvis.core.Atoms class several atomistic properties such as lattice parameters, density, packing fraction etc. can be calculated. The optimized lattice parameters generally compare well with DFT data except for the FFs where a particular phase was not trained during FF-fitting. Similarly, energetics in terms of convex hull plot is compared between DFT and FF results. An example is shown below:
 
Figure. Convex hull plot for Ni-Al system using a) DFT data, b) Force-field (Mishin Ni-Al potential). Filled blue circles show stable while the red squares show unstable structures. Here “fu” stands for formula unit.

Similarly, elastic tensor and derived properties were predicted using LAMMPS runs. 

 
Figure. Comparison of Voigt-bulk modulus data of force-fields and DFT data. The x-axis shows the bulk modulus for materials (for which FFs were available to us) from DFT calculations, while the y-axis shows the bulk modulus for corresponding materials using the force-field calculations.
Moreover, vacancy formation energies, surface energies can be compared between FF and DFT for FF quality assessment.
 
Fig. Vacancy and surface energies (up to a maximum Miller index of 3) for single elements obtained with all the FF available for such materials are compared to available DFT data.  
Phonon data from JARVIS-DFT and JARVIS-FF can be compared for a system to evaluate the phonon quality. However, it is important to note that in DFT there might be only conventional cell Gamma-point phonon data available whereas in JARVIS-FF we use supercell finite difference-based approach for obtaining phonon density of states and bandstructures. For a known stable material if the phonon bandstructure shows high negative values then it signifies the FF maynot be suitable to predict correct dynamical properties of the system.
