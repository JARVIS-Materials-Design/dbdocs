**JARVIS-DB Documentation**

This website contains information about the data available in the NIST-JARVIS ([https://jarvis.nist.gov/](https://jarvis.nist.gov/)) infrastructure. There are multiple apps and resources in JARVIS such as JARVIS-FF, JARVIS-DFT, JARVIS-ML, JARVIS-Tools, JARVIS-WannierTB, JARVIS-Heterostructure etc. Documentation on each of them contains the following sections:

- Overview,
- Statistics,
- Methodology,
- Property details, code example and post-processing,
- Data quality assessment,
- References.

**Outline**

1. **JARVIS-DFT**
  1. Overview
  2. Methodology
    1. Obtaining initial/input crystal structure
    2. Geometric and electronic optimization
    3. Functional selection
    4. K-point and cut-off convergence
    5. DFT+U
    6. Spin-orbit coupling
    7. Beyond DFT methods
    8. Density functional perturbation theory
    9. Finite-difference method
    10. Linear Optics
    11. Wannierization
    12. Boltzmann transport
    13. Partial charge density
  3. Property Details

1.3.1Atomic structure and derived properties

1.3.2 Formation energy

1.3.3 Exfoliation energies

1.3.4 Defect formation energies

1.3.5 Electronic density of density of states

1.3.6 Electronic Bandstructure

1.3.8 Frequency dependent dielectric function and optoelectronic properties

1.3.9 Static dielectric tensor

1.3.10 Piezoelectric tensor

1.3.11 Infrared intensity

1.3.12 Elastic tensor

1.3.13 Thermoelectric properties

1.3.14 Wannier tight binding Hamiltonians

1.3.15 Scanning tunneling microscopy images

    1. Electric field gradients
    2. FT convergence parameters

  1. Data Quality assessment table
  2. Statistics
  3. References

1. **JARVIS-FF**

2.1 Overview

2.2 Methodology

2.3 Property details and assessment

2.4 References

1. **JARVIS-ML**

3.1 Overview

3.2 Methodology

3.3 Quality assessment

3.4 References

**1 JARVIS-DFT**

**1.1 Overview**

JARVIS-DFT is a density functional theory-based database for ~40000 3D, ~1000 2D materials and around a million calculated properties. JARVIS-DFT mainly uses vdW-DF-OptB88 functional for geometry optimization. It also uses beyond-GGA approaches, including Tran-Blaha modified Becke-Johnson (TBmBJ) meta-GGA, PBE0, HSE06, DMFT, G0W0 for analyzing selective cases. In addition to hosting conventional properties such as formation energies, bandgaps, elastic constants, piezoelectric constants, dielectric constants, and magnetic moments, it also contains unique datasets, such as exfoliation energies for van der Waals bonded materials, spin-orbit coupling spillage, improved meta-GGA bandgaps, frequency-dependent dielectric function, spin-orbit spillage, spectroscopy limited maximum efficiency (SLME), infrared (IR) intensities, electric field gradient (EFG), heterojunction classifications, and Wannier tight-binding Hamiltonians. These datasets are compared to experimental results wherever possible, to evaluate their accuracy as predictive tools. JARVIS-DFT introduces protocols such as automatic k-point convergence that can be critical for obtaining precise and accurate calculation results.

![](RackMultipart20200915-4-htfc92_html_2aa90493fbdc4b54.png)

**1.2 Methodology**

**1.2.1 Obtaining initial/input crystal structure**

The initial/input crystal structures were obtained from the Inorganic Crystal Structure Database, Materials Project, OQMD, AFLOW, COD databases. These input structures are then subjected to a set of geometric and electronic optimizations (discussed below) following JARVIS-DFT workflow protocols. After the geometric optimization, several materials properties are calculated (discussed below). A typical flow chart is shown below:

![](RackMultipart20200915-4-htfc92_html_d3477193a83cec0f.png)

**1.2.2 Geometric and electronic optimization**

DFT calculations were carried out using the Vienna Ab-initio simulation package (VASP) software using the workflow given on our JARVIS-Tools github page ([https://github.com/usnistgov/jarvis](https://github.com/usnistgov/jarvis)). We use the projected augmented wave method .  Both the internal atomic positions and the lattice constants are allowed to relax in spin-unrestricted calculations until the maximal residual Hellmann–Feynman forces on atoms are smaller than 0.001 eV Å−1 and energy-tolerance of 10-7 eV with accurate precision setting (PREC=Accurate). Note that force convergence is very important for properties such as phonons, elastic constants etc. We use spin-polarized set-up with initial magnetic moment of default value 1 muB during geometric optimization. Also use conjugate gradient algorithm for ionic relaxation (IBRION=2).

**1.2.3 Functional selection**

We mainly used OptB88vdW method for our calculations, but we also carryout local density approximation (LDA) and generalized gradient approximation with Perdew-Burke-Ernzerhof (GGA-PBE)-based calculations for benchmarking purposes. OptB88vdW functional, has been shown to gives accurate lattice parameters for both van der Waals (vdW) and non-vdW solids. The exchange-correlation energy within vdW-DF is given by:

where is the exchange energy within the generalized gradient approximation (GGA) and is the correlation energy within the local-density approximation (LDA). The nonlocal correlation energy is given by:

where is a kernel function, and . The q0 is a function of n(r) and |∇n(r)|, and it is proportional to the gradient corrected LDA exchange-correlation energy per electron. This function controls the behavior of in the slowly varying as well as nonuniform density regions. It is noted that the use of the LDA correlation is motivated by the fact that vanishes in the uniform electron gas limit, and to avoid the possible double counting of the gradient correction contained in . Hence, the vdW-DF-optB88 is an example of the truly nonlocal-correlation functionals in the vdW-DF-method for approximating the vdW forces in regular DFT.

TBmBJ potential is used as a meta-GGA method for better predicting the bandgaps and dielectric function of a material. The TBmBJ potential is given by:

where _c_ is a system-dependent parameter, with _c_ = 1 corresponding to the Becke-Roussel (BR) potential, which was originally proposed to mimic the Slater potential, the Coulomb potential corresponding to the exact exchange hole. For bulk crystalline materials, Tran and Blaha proposed to determine _c_ by the following empirical relation:

With , Å1/2 and Vcell is the volume of the unit cell. The _c_-parameter was automatically determined in VASP through a self-consistent run.

**1.2.4 K-point and cut-off convergence**

We use the Monkhorst-Pack scheme to generate k-points, but after the generation, the grid is shifted so that one of the k-points lies on the Г-point. Here a &quot;k-points line density&quot; (L) is defined, which is related to the reciprocal lattice vectors by:

, , ,

 Where are the reciprocal lattice vectors, and   their norms.

We included the gamma-point because we were interested in computing quantities that require gamma-point contribution, such as optical transition for our optoelectronic database, gamma-point phonons for our elastic properties, finding multiferroic materials which have negative phonons at the gamma-point. The k-points are continuously stored in memory, to check that each of the new k-points generated by equation is unique. The k-points line density starts from length 0, with Г-point being the only k-point and is increased by 5 Å at each successive iteration if the difference between the energy computed with the new k-points and the one computed with previous k-points is more than the tolerance. After the convergence with a particular tolerance, we compute five extra points to further ensure the convergence. This procedure is repeated until convergence is found for all 5 extra points. A similar convergence procedure is carried out for the plane wave cut-off until the energy difference between successive iterations is less than the tolerance for 6 successive points. The plane wave cut-off is increased by 50 eV each time, starting from 500 eV. In both convergence procedures, we perform only single step electronic relaxation, i.e. no ionic relaxation is considered. When starting the cut-off energy convergence, we used a minimal k-point length of 10 Å. Similarly, for the k-point convergence we started with a cut-off of 500 eV. Note that complete ionic and electronic relaxation for determining converged parameters might be needed for very sophisticated calculation (such as Raman intensity calculation), but those calculations are beyond the scope of this work.

The plane wave cutoff is converged using the same procedure as for the 3D bulk materials. In this work, we used Gaussian smearing (with 0.01 eV parameter) which is recommended by several DFT codes, because it is less sensitive than other methods to determine partial occupancies for each orbital. This leads to an easier DFT-SCF convergence, especially when the materials are not apriori known to be a metal or insulator, which is always the case in this work. However, it is to be emphasized that, in principle, k-points and smearing parameters should be converged together, but this requires a very computationally expensive workflow. For this reason, we choose to converge k-points and cut-off only.

In addition to the above convergence procedure, we further increase the cut-off by 30% for elastic, dielectric and electric field gradient, piezoelectric tensors as the volume and shape of the material may change during the distortions.

**1.2.5 DFT+U**

DFT+U corrections are used only for special cases because it is somewhat tricky to get formation energies from the total energies obtained from such calculations. DFT+U calculations are used during magnetic topologic insulator search cases mainly. We generally use U=3.0 eV or complete U-scan (0 to 3 eV) for enumerating the effects of U-parameter. Refer to: [https://www.nature.com/articles/s41524-020-0319-4](https://www.nature.com/articles/s41524-020-0319-4)

**1.2.6 Spin-orbit coupling**

Spin-orbit coupling usually split states that are degenerate in a nonrelativistic description. We consider spin-orbit coupling (LSORBIT = .TRUE.) only during spin-orbit spillage for topological materials calculations but not during geometric optimization or any other major property calculations. We take into account only scalar relativistic effects. Accurate lattice-constants of materials with vdW bonding can play a critical role in predicting the correct topology, as emphasized. We calculate the spillage using the Perdew, Burke and Ernzerhof (PBE) functional and 1000/atom k-points, as well as by analyzing the spillage along a high-symmetry Brillouin-zone path, with a 600 eV plane-wave cut-off. We calculate the spin-orbit spillage, given by the following equation:

where

for wavefunctions without SOC and with SOC calculations.

**1.2.7 Beyond DFT methods**

JARVIS-DFT contains data using conventional local and semi-local methods, as well as beyond conventional DFT based methods such as meta-GGA (TBmBJ), G0W0, HSE06, PBE0, DMFT etc. BDFT methods are used to better predict electronic bandgaps, dielectric functions hence solar-cell efficiencies, as well magnetic moment (using DMFT) of a material. HSE06 and PBE0 were used to predict accurate bandgaps for exfoliable bulk and corresponding monolayer materials. We utilized two hybrid functionals: PBE0 and HSE06. In PBE0 the exchange energy is given by a 3:1 ratio mix of the PBE and Hartree–Fock exchange energies, respectively, while the correlation is completely given by the PBE correlation energy. In HSE (Heyd–Scuseria–Ernzerhof)32, the exchange is given by a screened Coulomb potential, to improve computational efficiency. An adjustable parameter (ω) controls how short range the interaction is. HSE06 is characterized by ω=0.2, while for a choice of ω=0 HSE becomes PBE0.

**1.2.8 Density functional perturbation theory**

We carry out the Density functional perturbation theory (DFPT) (as implemented in the VASP code) calculation on the standard conventional cell for each material. We determine the Born-effective charge, piezoelectric and dielectric (ionic+electronic parts) tensors and the phonon eigenvectors. DFPT is also used to calculate Infrared and Raman intensities for selected materials. As mentioned earlier, it is important to converge K-points, cut-off and parameters to get reliable results which is taken into account in the JARVIS-DFT. Given an insulating system with _N_ atoms per cell, with cell volume Ω0, atomic displacements um (m=(1…3N)), homogenous strain ηj (j=(1…6)), homogenous electric fields α (α=(x,y,z)), energy _E_, the force-constant matrix (), internal strain tensor (), the dielectric susceptibility (), Born dynamical effective charge (), piezoelectric stress tensor () and piezoelectric strain tensor () are calculated as follows (SI units are used throughout)65-69:

The dielectric constant can be derived from the dielectric susceptibility using:

The first term represents the electronic contribution and the second term the ionic contribution for DL and PZ constants respectively.

The PZ is a 3x6 tensor, the DL 3x3 and the BEC _N_x3x3 tensor. The IR intensity of phonon modes is calculated using:

where is the normalized vibrational eigenvector of the _n_th phonon mode of the _s_th atom in the unit cell, and α, β are the cartesian coordinates. is the Born effective charge tensor of _s_th atom (here we explicitly write both the cartesian indices of Z). These approaches are universal and have been already applied to various material classes. More details about the DFPT formalism can be found in elsewhere65,66.

**1.2.9 Finite-difference method**

The elastic tensor is determined by performing six finite distortions of the lattice and deriving the elastic constants from the strain-stress relationship. A set of strains **ε** where ε1, ε2, and ε3 are the normal strains and the others are the shear strains imposed on a crystal with lattice vectors **R** specified in Cartesian coordinates,

**R**

where _a __1_ is the x-component of the lattice vector , _b__ 2_ the y-component of the lattice vector , and so on. Corresponding to a set of strains discussed above, a set of stresses **σ** are determined with VASP code. The stress-strain can then be related by general Hooke&#39;s law

**σ = Cε**

where **C** is a 6x6 elastic constant matrix, which can be obtained by matrix-inverse operations. Finite difference method-based phonons were also obtained at Gamma-points only.

**1.2.10 Linear Optics**

To obtain the optical properties of the materials, we calculated the imaginary part of the dielectric function from the Bloch wavefunctions and eigenvalues (neglecting local field effects). We introduced three times as many empty conduction bands as valance bands. This treatment is necessary to facilitate proper electronic transitions. We choose 5000 energy grid points to have a sufficiently high resolution in dielectric function spectra. The imaginary part is calculated as:

where _e_ is electron charge, is the cell volume, _E_ the energy, ![](RackMultipart20200915-4-htfc92_html_5b6874a2aedb81f5.jpg)is the Fermi-weight of each _k_-point, ![](RackMultipart20200915-4-htfc92_html_6a50251f1cecb5e4.jpg)are unit vectors along the three Cartesian directions, is the cell-periodic part of the pseudo-wavefunction for band _n_ and k-point , _q_ stands for the Bloch vector of an incident wave, _c_ and _v_ stand for conduction and valence bands, stands for eigenvalues of the corresponding bands respectively. The matrix elements on the right side of Eq. (1) capture the transitions allowed by symmetry and selection rules. The real part of the dielectric tensor  is obtained by the usual Kramers-Kronig transformation:

where _P_ denotes the principle value, and _η_ is the complex shift parameter taken as 0.1. Moreover, as the dielectric function is a tensorial quantity, we use the crystallographic average of the dielectric function (written as and ), obtained by diagonalizing the dielectric tensor for each energy and averaging the diagonal elements.

Using, and the absorption coefficient is defined as:

where _c_ is the speed of light.

Next, the SLME ( is defined as the ratio of the maximum output power density (_P__max_) and the total incident solar energy density (_P__in_). _P__max_ is obtained by numerically maximizing the product of current density _J_ and voltage _V_.

Assuming the solar cell at temperature _T_ behaves as an ideal diode and is illuminated under the photon flux _I__sun_, _J_ and _V_ follow the following equation:

where _e_ is the elementary charge, _V_ the potential over the absorber layer and _k_ is Boltzmann&#39;s constant. The first term is the short-circuit current density _J__sc_ given by:

where and _a(E)_ is the photon absorptivity, is the AM1.5G solar spectrum. The _a(E)_ depends on the absorption coefficient () (Eq. (3)) and thickness (_L_) of the material.

The coefficient of the second term in Eq. (5) is the reverse saturation current (_J__0_), which corresponds to the total (radiative and non-radiative) electron-hole recombination current at equilibrium in the dark:

Here, _f __r_ is defined as the fraction of the radiative recombination current. For the SLME, _f__ r_ is approximated using:

Where is the fundamental and is the direct allowed bandgap of a material.

Following the principle of detailed balance, the rates of emission and absorption through cell surfaces must be equal in equilibrium in the dark. Hence, can be calculated from the rate at which black-body photons from the surrounding thermal bath are absorbed through the front surface, given by:

where _I__bb_ is the black-body spectrum at temperature _T._ Both the solar spectrum and black-body spectrum are expressed in terms of the photon flux.

In order to maximize the power density, Eq. can be re-written as:

Therefore, the material-property related inputs in calculating the SLME are and T. In this work, we assume material thickness (_L_) as 500 nm and temperature (_T_) as 300 K.

VASP uses a complex shift (CSHIFT) in the Kramers-Kronig relation to calculate the real part of the dielectric tensor, and also determines the corresponding imaginary part for consistency.

**1.2.11 Wannierization**

We use Wannier90 to construct Maximally-Localized Wannier Functions (MLWF) based TB-Hamiltonians. The basic formalism of Wannierization is well-established. We briefly review some aspects here, interested readers can see longer discussions in 5,16. For a set of Bloch eigenvectors , a general set of WFs (n=1...N) can be written as:

wherelabels the unit cell of the WF, _V_ is the volume of the unit cell, and is an arbitrary unitary matrix. To construct maximally-localized WFs, is chosen to minimize the following spread functional:

where and . The minimization proceeds iteratively, based on an initial guess of localized orbitals.

For the case of interest in this work, where we wish to describe both the valence and conduction bands near the Fermi level, it is necessary to first select a set of bands to Wannierize, which includes separating the conduction bands from the free-electron-like bands that generally overlap with them in energy. The procedure to determine this localized subspace of Bloch wavefunctions proceeds similarly to minimization described above, where after an initial guess, the subspace is iteratively updated in order to minimize the spread function. After this initial disentanglement step, the Wannierization of the selected subspace proceeds as described above.

Due to the iterative non-linear minimization employed during both the disentanglement and Wannierization steps, the localization and utility of the final Wannier functions depend in practice on the initial choice of orbitals that are used to begin the disentanglement procedure, and which are then used as the initial guess for the Wannierization. Our initial guesses consist of a set of atomic orbitals we have chosen to describe all the chemically relevant orbitals for each element in typical elemental systems and compounds. We provide the list of the orbitals we select for each element in Table.

_Table: Semi-core states used during Wannierization. Note that several lanthanides and actinides are avoided._

| **Element** | **POTCARs** | **N\_electrons** | **Excluded\_electrons** | **Projections** | **N\_Wan\_projections** |
| --- | --- | --- | --- | --- | --- |
| **Ag** | Ag | 11 | 0 | s\_d | 6 |
| **Al** | Al | 3 | 0 | s\_p | 4 |
| **Ar** | Ar | 8 | 0 | s\_p | 4 |
| **As** | As | 5 | 0 | s\_p | 4 |
| **Au** | Au | 11 | 0 | s\_d | 6 |
| **B** | B | 3 | 0 | s\_p | 4 |
| **Ba** | Ba\_sv | 2 | 8 | s\_d | 6 |
| **Be** | Be\_sv | 2 | 2 | s\_p | 4 |
| **Bi** | Bi | 5 | 0 | s\_p | 4 |
| **Br** | Br | 7 | 0 | s\_p | 4 |
| **C** | C | 4 | 0 | s\_p | 4 |
| **Ca** | Ca\_sv | 2 | 8 | s\_d | 6 |
| **Cd** | Cd | 2 | 10 | s\_d | 6 |
| **Ce** | Ce | 4 | 8 | f\_d\_s | 13 |
| **Cl** | Cl | 5 | 2 | p | 3 |
| **Co** | Co | 9 | 0 | s\_d | 6 |
| **Cr** | Cr\_pv | 6 | 6 | s\_d | 6 |
| **Cs** | Cs\_sv | 1 | 8 | s\_d | 6 |
| **Cu** | Cu\_pv | 11 | 6 | s\_d | 6 |
| **Dy** | Dy\_3 | 12 | 3 | s\_f | 8 |
| **Er** | Er\_3 | 14 | 5 | f\_s | 8 |
| **Eu** | Eu | 9 | 8 | f\_s | 8 |
| **F** | F | 5 | 2 | p | 3 |
| **Fe** | Fe\_pv | 8 | 6 | s\_d | 6 |
| **Ga** | Ga\_d | 3 | 10 | s\_p | 4 |
| **Gd** | Gd | 10 | 8 | f\_d\_s | 8 |
| **Ge** | Ge\_d | 4 | 10 | s\_p | 4 |
| **H** | H | 1 | 0 | s | 1 |
| **He** | He | 2 | 0 | s | 1 |
| **Hf** | Hf\_pv | 4 | 6 | s\_d | 6 |
| **Hg** | Hg | 12 | 0 | s\_p\_d | 9 |
| **I** | I | 7 | 0 | s\_p | 4 |
| **In** | In\_d | 3 | 10 | s\_p | 4 |
| **Ir** | Ir | 9 | 0 | s\_d | 6 |
| **K** | K\_sv | 1 | 8 | s\_d | 6 |
| **Kr** | Kr | 8 | 0 | s\_p | 4 |
| **La** | La | 3 | 8 | s\_d\_f | 13 |
| **Li** | Li\_sv | 1 | 2 | s | 1 |
| **Lu** | Lu\_3 | 17 | 8 | f\_d\_s | 13 |
| **Mg** | Mg\_pv | 2 | 6 | s\_p | 4 |
| **Mn** | Mn\_pv | 7 | 6 | s\_d | 6 |
| **Mo** | Mo\_pv | 6 | 6 | s\_d | 6 |
| **N** | N | 3 | 2 | p | 3 |
| **Na** | Na\_pv | 1 | 6 | s\_p | 4 |
| **Nb** | Nb\_pv | 5 | 6 | s\_d | 6 |
| **Nd** | Nd\_3 | 6 | 5 | f\_s | 8 |
| **Ne** | Ne | 8 | 0 | s\_p | 4 |
| **Ni** | Ni\_pv | 10 | 6 | s\_d | 6 |
| **O** | O | 4 | 2 | p | 3 |
| **Os** | Os\_pv | 8 | 6 | s\_d | 6 |
| **P** | P | 5 | 0 | s\_p | 4 |
| **Pb** | Pb\_d | 4 | 10 | s\_p | 4 |
| **Pd** | Pd | 10 | 0 | s\_d | 6 |
| **Pt** | Pt | 10 | 0 | s\_d | 6 |
| **Rb** | Rb\_sv | 1 | 8 | s\_d | 6 |
| **Re** | Re\_pv | 7 | 6 | s\_d | 6 |
| **Rh** | Rh\_pv | 9 | 6 | s\_d | 6 |
| **Ru** | Ru\_pv | 8 | 6 | s\_d | 6 |
| **S** | S | 4 | 2 | p | 3 |
| **Sb** | Sb | 5 | 0 | s\_p | 4 |
| **Sc** | Sc\_sv | 3 | 8 | s\_d | 6 |
| **Se** | Se | 6 | 0 | s\_p | 4 |
| **Si** | Si | 4 | 0 | s\_p | 4 |
| **Sm** | Sm\_3 | 8 | 3 | f\_s | 8 |
| **Sn** | Sn\_d | 4 | 10 | s\_p | 4 |
| **Sr** | Sr\_sv | 2 | 8 | s\_d | 6 |
| **Ta** | Ta\_pv | 5 | 6 | s\_d | 6 |
| **Tb** | Tb\_3 | 9 | 0 | f\_s | 8 |
| **Tc** | Tc\_pv | 7 | 6 | s\_d | 6 |
| **Te** | Te | 6 | 0 | s\_p | 4 |
| **Th** | Th | 4 | 8 | d\_s | 6 |
| **Ti** | Ti\_pv | 4 | 6 | s\_d | 6 |
| **Tl** | Tl\_d | 3 | 10 | s\_p | 4 |
| **U** | U | 9 | 5 | f\_s | 8 |
| **V** | V\_pv | 5 | 6 | s\_d | 6 |
| **W** | W\_pv | 6 | 0 | s\_d | 6 |
| **Xe** | Xe | 8 | 0 | s\_p | 4 |
| **Y** | Y\_sv | 3 | 8 | s\_d | 6 |
| **Zn** | Zn | 12 | 0 | s\_p\_d | 9 |
| **Zr** | Zr\_sv | 4 | 8 | s\_d | 6 |

For many specific materials, it may be possible to select a smaller set of orbitals while still maintaining high-quality WFs that describe the bands of interest; however, our fairly inclusive set of orbitals is able Wannierize nearly all compounds in a high-throughput manner without human intervention. Because most applications of WFs are computationally inexpensive compared to the DFT calculations used to construct the WFs, in practice, our larger Wannier basis has only minimal computational cost. However, it is necessary to have enough empty bands in the underlying DFT calculation such that any empty orbitals chosen are included in the Bloch basis. We do not include any semicore orbitals in our Wannier basis, as they are generally well-separated in energy from the valence orbitals and are not necessary to describe bands near the Fermi level, and we exclude semicore bands from the disentanglement.

During the disentanglement step, it is possible to choose an energy range that is included exactly (&quot;the frozen window&quot;), with the result that the Wannier band structure will exactly match the DFT band structure in this energy range and at the grid of k-points used in the Wannierization. We use a frozen window of 2 eV around the Fermi-energy. This window ensures that bands near the Fermi level are well described by the WTBH. For cases where the original WFs were unsatisfactory (see below), we found that lowering the lower bound of this window to include all of the valence bands often improves that WTBH, which we use as a second Wannierization setting. In order to validate our WTBH, we calculate the maximum absolute difference () between the Wannier and DFT eigenvalues within an energy range of ± 2eV around the Fermi level:

As discussed above, at the grid of k-points used in the construction of the WFs and within the frozen window, the eigenvalues should match exactly by construction. Therefore, to test the WTBH, we evaluate Eq. on the dense lines of k-points used to generate our band structure plots, which tests the WFs out of sample. A weakness of this evaluation method is that highly dispersive energy bands (high ) can result in high values even if the WTBH is of good quality because any slight shift in the _k_-direction of a dispersive band will result in a large energy error. We consider that systems with less than 0.1 eV to useful for most applications, and we provide data for the user to evaluate individual WTBH for their own applications.

**1.2.12 Boltzmann transport**

The transport properties were calculated using the Boltzmann transport equation (BTE) implemented in the BoltzTrap code. The BTE is used to investigate the non-equilibrium behavior of electrons and holes by statistically averaging all possible quantum states given by the equation below:

.

_f_ is the electron distribution, which is a function of state _k_, temperature _T_ and time _t_, and _v(k)_ are the electron group velocities. The three terms on the right-hand side of Eq. refer, respectively, to the temporal rate of change of _f_ due to all scattering processes, the rate of change of _f_ due to external forces, and the diffusion from the carrier density gradient. If the external forces consist only of a low electric field, E, and no magnetic field, B, such that then the low-filed BTE is given by:

.

Now, _f_ can be described as a first-order (linear) perturbation from the (equilibrium) Fermi-Dirac distribution, _f__0_, due to scattering

where

where the dependence of on _k_ is given by the electronic band structure, and the various

scattering terms and time dependence are lumped into the electronic relaxation time,. The computation of the relaxation time is very computationally expensive, especially in a high-throughput context.

**1.2.13 Partial charge density**

The partial charge densities are used to calculate scanning tunneling microscopy images. The surface charge and probability densities are calculated by integrating the local density of states function (ILDOS) over an energy range of ±0.5 eV from the conduction band minima (CBM) to Fermi energy (EF) and valence band maxima (VBM) to Fermi energy (EF). The STM images are calculated using Tersoff-Hamann approach:

where is the tunneling current, the eigenvector and the eigenenergy due to the state μ and is the Fermi-energy.

  1. **Property Details**

**1.3.1 Atomic structure and derived properties**

Atomic structure is generally defined using lattice parameter and atomic positions (xyz coordinate) and element information. Experimentally, this information is obtained with experimental techniques such as X-ray diffraction, EBSD, Neutron diffraction, NMR etc. After geometric optimization following the JARVIS-DFT protocol several atomic properties are calculated. Some of these properties are lattice-constants/parameters (lengths and angles), density, packing fraction, spacegroup, crystal system, radial distribution, angular distribution, dihedral distribution, X-ray diffraction information. This can be done with the help of JARVIS-Tools jarvis.core.Atoms class. Example:

![](RackMultipart20200915-4-htfc92_html_49e270ebaec91b2.png)

![](RackMultipart20200915-4-htfc92_html_f279d1dfa772ee41.png)

![](RackMultipart20200915-4-htfc92_html_1fb7b3fa53e34903.png)

_Table Mean absolute error (MAE) and root-mean-squared error (RMSE) in a, b and c crystallographic directions computed for all materials in our database with respect to experimental data (ICSD data). To facilitate comparison between the functionals, both MAE and RMSE have been computed for all materials, only for predicted vdW bonded materials and only for predicted non-vdW bonded materials, using Material&#39;s project PBE and JARVIS-DFT OPT functional._

|
 | **#Mats.** | **MAE (_a_)** | **MAE (_b_)** | **MAE (_c_)** | **RMSE (_a_)** | **RMSE (_b_)** | **RMSE (_c_)** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **OPT (All)** | 10052 | 0.11 | 0.11 | 0.18 | 0.29 | 0.30 | 0.58 |
| **PBE (All)** | 10052 | 0.13 | 0.14 | 0.23 | 0.30 | 0.29 | 0.61 |
| **OPT (vdW)** | 2241 | 0.20 | 0.21 | 0.44 | 0.44 | 0.44 | 0.99 |
| **PBE (vdW)** | 2241 | 0.26 | 0.29 | 0.62 | 0.45 | 0.51 | 1.09 |
| **OPT (non-vdW)** | 7811 | 0.08 | 0.08 | 0.11 | 0.23 | 0.24 | 0.39 |
| **PBE (non-vdW)** | 7811 | 0.09 | 0.09 | 0.12 | 0.22 | 0.25 | 0.36 |

**1.3.2 Formation energy**

The enthalpy of formation is the standard reaction enthalpy for the formation of the compound from its elements (atoms or molecules) in their most stable reference states. We report formation energies at 0K only based on energies directly obtained from the DFT especially with OptB88vdW functional. We have calculated the elemental energies (in respective crystalline forms) and treat those energies as the chemical potential. These chemical potentials are subtracted from the total energies of the system to predict the formation energies. If functional other than OptB88vdW is used then the unary energy needs to be provided.

![](RackMultipart20200915-4-htfc92_html_4187f13b2c8ec9f7.png)

**1.3.3 Exfoliation energies**

The exfoliation energy for 2D materials is computed as equation:

Here, _E __1L_ and _E__ bulk_are the energies of the single layer and 3D bulk materials and _N __1L_ and _N__ bulk_are the number of the atoms in the single layer and bulk systems respectively.

![](RackMultipart20200915-4-htfc92_html_22647b4c2c9386bb.png)

_Fig. Exfoliation energy distribution for the materials based on the energy difference of bulk and layered materials phase. Most of the materials have E__f_ _\&lt; 200 meV/atom suggesting experimental fabrication of layered structures. Materials inside each bins are grouped and discussed below. Most of the materials are found in the 60-100 meV range._

**1.3.4 Defect formation energies**

Defect formation energies are available for a few 2D materials. For neutral defect formation energies, we generate a at least 11 Angs. Cell for 2D materials and create point defects using unique Wyckoff information. The total energy of the defect system, the energy of the perfect 2D crystal and the elemental chemical potential is then used to predict the defect formation energies. Currently this property is available for OptB88vdW functional only.

![](RackMultipart20200915-4-htfc92_html_9d4043e3ea288b2.png)

**1.3.5 Electronic density of density of states**

EDOS is a spectrum of the number of allowed electron energy levels (occupied and unoccupied) against energy interval (generally in electron volt). EDOS is generally calculated based on a dense k-point grid with a smoothening/smearing setup. EDOS is generally rescaled with respect to the Fermi-energy of the system. A high value for the DOS represents a high number for the electronic states that can be occupied. Experimental methods such as scanning tunneling microscopy can be used obtain EDOS. One of the obvious properties that can be calculated from EDOS is the bandgap. EDOS can be based on total density of states or further resolved into atom-projected and element-projected density of states. The projected EDOS provides information about the particular electronic orbitals (say s,p,d,f) and specific atom (say Si) contributing towards a particular energy level. We calculate EDOS for both semi-local, HSE06 and PBE0 methods. At least, DOS is provided at local or semilocal materials for all the materials. An example is shown below:

![](RackMultipart20200915-4-htfc92_html_be7fd68ec9622ac8.png)

An example of calculating EDOS using JARVIS-tools is given below:

![](RackMultipart20200915-4-htfc92_html_dc21e121fd8f9ad9.png)

**1.3.6 Electronic Bandstructure**

Following the laws of quantum mechanics only certain energy levels are allowed and others could be forbidden. These levels can be discrete or spilt called bands. Electronic bands depend on crystal momentum leading to bandstructure. The allowed states which are filled (upto the Fermi level) are valence bands while the unoccupied bands are conduction bands. The difference in energy between the conduction and valence band gives rise to the bandgaps. The materials with no bandgap are termed metals, with low bandgap semiconductors and with high bandgaps insulators. We calculate bandstructure using local, semi-local as well as hybrid functionals. An example is shown below:

![](RackMultipart20200915-4-htfc92_html_9c73a0d042d73191.png)

band\_info=vrun.get\_bandstructure(kpoints\_file\_path=&#39;KPOINTS&#39;)

A comparison of bandgaps for selected materials for which we could find experimental results is shown below.

_Table: Comparison of bandgaps obtained from OPT functional and MBJ potential schemes compared with experimental results and DFT data available in different databases. Materials, space-group (SG), Inorganic Crystal Structure Database (ICSD#) id, Materials-Project (MP#) id, JARVIS-DFT id (JV#), bandgap from MP (MP), bandgap from AFLOW, bandgap from OQMD, our OptB88vdW bandgap (OPT), Tran-Blah modified Becke-Johnson potential bandgap (MBJ), Heyd-Scuseria-Ernzerhof (HSE06) and experimental bandgaps (eV) data are shown. Experimental data were obtained from_ _18,21,46,47__. MAE denotes the mean absolute error, while SC is the Spearman&#39;s coefficient._

| **Mats.** | **SG** | **ICSD#** | **MP#** | **JV#** | **MP** | **AFLOW** | **OQMD** | **OPT** | **MBJ** | **HSE06** | **Exp.** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C | Fd-3m | 28857 | 66 | 91 | 4.11 | 4.12 | 4.4 | 4.46 | 5.04 | 5.26 | 5.5 |
| Si | Fd-3m | 29287 | 149 | 1002 | 0.61 | 0.61 | 0.8 | 0.73 | 1.28 | 1.22 | 1.17 |
| Ge | Fd-3m | 41980 | 32 | 890 | 0.0 | 0 | 0.4 | 0.01 | 0.61 | 0.82 | 0.74 |
| BN | P63/mmc | 167799 | 984 | 17 | 4.48 | 4.51 | 4.4 | 4.46 | 6.11 | 5.5 | 6.2 |
| AlN | P63mc | 31169 | 661 | 39 | 4.06 | 4.06 | 4.5 | 4.47 | 5.20 | 5.49 | 6.19 |
| AlN | F-43m | 67780 | 1700 | 7844 | 3.31 | 3.31 | - | 3.55 | 4.80 | 4.55 | 4.9 |
| GaN | P63mc | 34476 | 804 | 30 | 1.74 | 1.91 | 2.1 | 1.94 | 3.08 | 3.15 | 3.5 |
| GaN | F-43m | 157511 | 830 | 8169 | 1.57 | 1.75 | - | 1.79 | 2.9 | 2.85 | 3.28 |
| InN | P63mc | 162684 | 22205 | 1180 | 0.0 | 0.0 | - | 0.23 | 0.76 | - | 0.72 |
| BP | F-43m | 29050 | 1479 | 1312 | 1.24 | 1.25 | 1.4 | 1.51 | 1.91 | 1.98 | 2.1 |
| GaP | F-43m | 41676 | 2490 | 1393 | 1.59 | 1.64 | 1.7 | 1.48 | 2.37 | 2.28 | 2.35 |
| AlP | F-43m | 24490 | 1550 | 1327 | 1.63 | 1.63 | 1.7 | 1.79 | 2.56 | 2.30 | 2.50 |
| InP | F-43m | 41443 | 20351 | 1183 | 0.47 | 0.58 | 0.7 | 0.89 | 1.39 | 1.43 | 1.42 |
| Mats. | SG | ICSD# | MP# | JV# | MP | AFLOW | OQMD | OPT | MBJ | HSE | Exp. |
| AlSb | F-43m | 24804 | 2624 | 1408 | 1.23 | 1.23 | 1.4 | 1.32 | 1.77 | 1.80 | 1.69 |
| InSb | F-43m | 24519 | 20012 | 1189 | 0.0 | 0.0 | 0.0 | 0.02 | 0.80 | 0.45 | 0.24 |
| GaAs | F-43m | 41674 | 2534 | 1174 | 0.19 | 0.30 | 0.8 | 0.75 | 1.32 | 1.40 | 1.52 |
| InAs | F-43m | 24518 | 20305 | 97 | 0.0 | 0.0 | 0.3 | 0.15 | 0.40 | 0.45 | 0.42 |
| BAs | F-43m | 43871 | 10044 | 7630 | 1.2 | 1.2 | 1.4 | 1.42 | 1.93 | 1.86 | 1.46 |
| **MoS2** | P63/mmc | 24000 | 2815 | 54 | 1.23 | 1.25 | 1.3 | 0.92 | 1.33 | 1.49 | 1.29 |
| **MoSe2** | P63/mmc | 49800 | 1634 | 57 | 1.42 | 1.03 | 1.0 | 0.91 | 1.32 | 1.40 | 1.11 |
| **WS2** | P63/mmc | 56014 | 224 | 72 | 1.56 | 1.29 | 1.4 | 0.72 | 1.51 | 1.6 | 1.38 |
| **WSe2** | P63/mmc | 40752 | 1821 | 75 | 1.45 | 1.22 | 1.2 | 1.05 | 1.44 | 1.52 | 1.23 |
| **Al2O3** | R-3c | 600672 | 1143 | 32 | 5.85 | 5.85 | 6.3 | 6.43 | 7.57 | 8.34 | 8.8 |
| CdTe | F-43m | 31844 | 406 | 23 | 0.59 | 0.71 | 1.1 | 0.83 | 1.64 | 1.79 | 1.61 |
| SnTe | Fm-3m | 52489 | 1883 | 7860 | 0.04 | 0.25 | 0.3 | 0.04 | 0.16 | 0.17 | 0.36 |
| SnSe | Pnma | 60933 | 691 | 299 | 0.52 | - | 0.6 | 0.71 | 1.25 | 0.89 | 0.90 |
| MgO | Fm-3m | 9863 | 1265 | 116 | 4.45 | 4.47 | 5.3 | 5.13 | 6.80 | 7.13 | 7.83 |
| CaO | Fm-3m | 26959 | 2605 | 1405 | 3.63 | 3.64 | 3.8 | 3.74 | 5.29 | 5.35 | 7.0 |
| CdS | P6\_3mc | 31074 | 672 | 95 | 1.2 | 1.25 | 1.4 | 1.06 | 2.61 | - | 2.5 |
| CdS | F-43m | 29278 | 2469 | 8003 | 1.05 | 1.19 | 1.4 | 0.99 | 2.52 | 2.14 | 2.50 |
| CdSe | F-43m | 41528 | 2691 | 1192 | 0.51 | 0.64 | 1.0 | 0.79 | 1.84 | 1.52 | 1.85 |
| MgS | F-43m | 159401 | 1315 | 1300 | 2.76 | 3.39 | 3.6 | 2.95 | 4.26 | 4.66 | 4.78 |
| MgSe | Fm-3m | 53946 | 10760 | 7678 | 1.77 | 1.77 | 1.8 | 2.12 | 3.37 | 2.74 | 2.47 |
| Mats. | SG | ICSD# | MP# | JV# | MP | AFLOW | OQMD | OPT | MBJ | HSE | Exp. |
| MgTe | F-43m | 159402 | 13033 | 7762 | 2.32 | 2.32 | 2.5 | 2.49 | 3.49 | 3.39 | 3.60 |
| BaS | Fm-3m | 30240 | 1500 | 1315 | 2.15 | 2.15 | 2.4 | 2.15 | 3.23 | 3.11 | 3.88 |
| BaSe | Fm-3m | 43655 | 1253 | 1294 | 1.95 | 1.95 | 2.9 | 1.97 | 2.85 | 2.79 | 3.58 |
| BaTe | Fm-3m | 29152 | 1000 | 1267 | 1.59 | 1.59 | 1.7 | 1.61 | 2.15 | 2.31 | 3.08 |
| **TiO2** | P42/mnm | 9161 | 2657 | 5 | 1.78 | 2.26 | 1.8 | 1.77 | 2.07 | 3.34 | 3.30 |
| **TiO2** | I41/amd | 9852 | 390 | 104 | 2.05 | 2.53 | 2.0 | 2.02 | 2.47 | - | 3.4 |
| **Cu2O** | Pn-3m | 26183 | 361 | 1216 | 0.5 | - | 0.8 | 0.13 | 0.49 | 1.98 | 2.17 |
| **CuAlO2** | R-3m | 25593 | 3748 | 1453 | 1.8 | 2.0 | 2.4 | 2.06 | 2.06 | - | 3.0 |
| **ZrO2** | P21/c | 15983 | 2858 | 113 | 3.47 | 3.56 | 4.0 | 3.62 | 4.21 | - | 5.5 |
| **HfO2** | P21/c | 27313 | 352 | 9147 | 4.02 | 4.02 | 4.5 | 4.12 | 5.66 | - | 5.7 |
| CuCl | F-43m | 23988 | 22914 | 1201 | 0.56 | 1.28 | 0.8 | 0.45 | 1.59 | 2.37 | 3.4 |
| **SrTiO3** | Pm-3m | 23076 | 5229 | 8082 | 2.1 | 2.29 | 1.8 | 1.81 | 2.30 | - | 3.3 |
| ZnS | F-43m | 41985 | 10695 | 1702 | 2.02 | 2.67 | 2.4 | 2.09 | 3.59 | 3.30 | 3.84 |
| ZnSe | F-43m | 41527 | 1190 | 96 | 1.17 | 1.70 | 1.5 | 1.22 | 2.63 | 2.37 | 2.82 |
| ZnTe | F-43m | 104196 | 2176 | 1198 | 1.08 | 1.48 | 1.5 | 1.07 | 2.23 | 2.25 | 2.39 |
| SiC | F-43m | 28389 | 8062 | 8158 | 1.37 | 1.37 | 1.5 | 1.62 | 2.31 | - | 2.42 |
| LiF | Fm-3m | 41409 | 1138 | 1130 | 8.72 | 8.75 | 11.0 | 9.48 | 11.2 | - | 14.2 |
| KCl | Fm-3m | 18014 | 23193 | 1145 | 5.03 | 5.05 | 5.3 | 5.33 | 8.41 | 6.53 | 8.50 |
| AgCl | Fm-3m | 56538 | 22922 | 1954 | 0.95 | 1.97 | 1.1 | 0.93 | 2.88 | 2.41 | 3.25 |
| AgBr | Fm-3m | 52246 | 23231 | 8583 | 0.73 | 1.57 | 0.9 | 1.00 | 2.52 | 2.01 | 2.71 |
| AgI | Fm-3m | 52361 | 22919 | 8566 | 0.77 | 1.98 | 1.4 | 0.39 | 2.08 | 2.48 | 2.91 |
| **MAE1** | - | - | - | - | 1.45 | 1.23 | 1.14 | 1.33 | 0.51 | 0.41 | - |
| **MAE2** | - | - | - | - | 1.39 | 1.19 | 1.09 | 1.27 | 0.43 | 0.42 | - |
| S.C. | - | - | - | - | 0.81 | 0.94 | 0.88 | 0.84 | 0.94 | 0.97 | - |

1. MAE calculated with respect to experiment using all available data for each method
2. MAE calculated with respect to experiment using only data for materials that have results available in all three DFT methods.

**1.3.7 Spin-orbit spillage and topological properties**

Topological materials driven by spin-orbit coupling have different bandgaps with/without spin orbit coupling effects. Spillage is calculated by comparing the wavefunctions of these two bandstructures. We calculate the spillage values for a large set of low bandgap and high atomic mass materials, with a spillage value of 0.5 as a threshold to screen for potential topologically non-trivial materials. The threshold denotes number of band-inverted electrons. This criterion can be used for predicting non-trivial topological behavior of metals, semiconductors and both perfect and defective structures. Atypical example of SOC spillage and band-inversion for Bi2Se3 is shown below:

![](RackMultipart20200915-4-htfc92_html_d7b4f31c5b6b8908.png)

![](RackMultipart20200915-4-htfc92_html_688573fd436054eb.png)

Flow chart for finding 3D non-magnetic topological materials is given below:

![](RackMultipart20200915-4-htfc92_html_722e338f85e37043.png)

A typical flowchart for 2D topological materials is given below:

![](RackMultipart20200915-4-htfc92_html_3eceaa1a3321154e.png)

Periodic table trend for high spillage materials are shown below:

![](RackMultipart20200915-4-htfc92_html_3309e9143529594.png)

**1.3.8 Frequency dependent dielectric function and optoelectronic properties**

We use linear optics theory at semilocal and meta-GGA levels to calculate frequency dependent real and imaginary part of dielectric function. The meta-GGA based data should predict better bandgaps hence dielectric function. Note that ionic contributions are ignored in such calculations. Several important properties such as absorption coefficient, electron energy loss spectrum, reflectivity, solar cell efficiency etc. can be calculated from the frequency dependent dielectric function.

![](RackMultipart20200915-4-htfc92_html_e9bd80299c2bf9a3.png)

_Fig. Imaginary part of the dielectric function in the x-direction_ _obtained from OPT, MBJ_ _and experiments._ _a) 1T-SnSe__2_ _(), b) 2H-MoS__2_ _(P6 __3__ /mmc), c) Si (), d) Ge (), e) GaAs ( and f) InP ._

![](RackMultipart20200915-4-htfc92_html_a70c8bc441f4f15b.png)

![](RackMultipart20200915-4-htfc92_html_a65646247d466b0e.png)

Using the absorption coefficient and bandgap (indirect and direct gaps), we calculate theoretical solar cell efficiencies of a material.

![](RackMultipart20200915-4-htfc92_html_1cf14bd28b33e86a.png)

_Fig. Summary of SLME data. a) SLME-distribution of all the materials in the database, b) SLME (ɳ) vs TBmBJ fundamental bandgap for high-SLME materials, c) colormap of SLME values with the direct versus indirect bandgaps for high-SLME materials, d) SLME vs average effective mass of electrons, e) dimensionality distribution in terms of 3D-bulk, 2D-bulk, 1D-bulk and 0D-bulk materials in the database, f) number of species distribution for high SLME materials, g) space-group distributions (1-230) for high-SLME materials, h and i) compositional-prototype and crystal-system distributions for high-SLME materials._

![](RackMultipart20200915-4-htfc92_html_6ef7d43bcecf5bc9.png)

_Fig. Periodic table trends of high SLME materials. The elements in a material are weighed 1 or 0 if the material has high or low-SLME. Then the probability of finding the element in a high-SLME material is calculated._

_Table Bandgap and SLME properties of a selection of materials with TBmBJ and G __0__ W __0_ _methods in DFT to evaluate uncertainty in predictions. Here E__ g_ _denotes the bandgap in eV and ɳ the calculated SLME in percentage._

| **Materials** | **JID** | **E**** g **** (TBmBJ) **|** E ****g** **(G****0 ****W**** 0****)** | **E**** g **** (G ****0**** W ****0**** +SOC) **|** Ƞ (TBmBJ) **|** Ƞ (G ****0**** W ****0**** ) **|** Ƞ (G ****0**** W ****0**** +SOC)** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **CuBr** | 5176 | 1.9 | 2.01 | 2.09 | 25.6 | 22.74 | 21.04 |
| **AuI** | 3849 | 2.1 | 2.34 | 2.20 | 11.4 | 8.83 | 11.86 |
| **SiAs** | 4216 | 1.6 | 1.36 | 1.33 | 26.1 | 23.85 | 23.20 |
| **BiTeBr** | 8781 | 1.90 | 1.52 | 0.79 | 25.2 | 32.15 | 26.11 |
| **TlPt**** 2 ****S**** 3** | 4630 | 1.30 | 1.45 | 1.35 | 32.70 | 30.99 | - |
| **MAD** | - | - | 0.22 | 0.34 | - | 3.23 | 2.21 |

**1.3.9 Static dielectric tensor**

Dielectric materials are important components in many electronic devices such as capacitors, field-effect transistors computer memory (DRAM), sensors and communication circuits. Both high and low-value DL materials have applications in different technologies. Static dielectric constants with both ionic and electronic contributions are obtained using DFPT method. The electronic part of the dielectric constant is generally in close agreement with that obtained from the linear optics method mentioned above.

![](RackMultipart20200915-4-htfc92_html_d8c5bc355308523.png)

_Fig. Flow-chart portraying different steps for the DFPT calculations._

_Table. Comparison of static dielectric constant for DFPT, MBJ and experiment. Experimental data were obtained from. MBJ data were obtained from our optoelectronic property database._

| **Materials** | **JID** | **DFPT** | **MBJ** | **Experiment** |
| --- | --- | --- | --- | --- |
| **MoS2** | 54 | _ε_11 =15.56 | _ε_11=15.34 | _ε_11=17.0 |
| **MoSe2** | 57 | _ε_11=16.90 | _ε_11=16.53 | _ε_11=18.0 |
| **MoTe2** | 60 | _ε_11=21.72 | _ε_11=18.74 | _ε_11=20.0 |
| **WS2** | 72 | _ε_11=13.91 | _ε_11=13.95 | _ε_11=11.5 |
| **WSe2** | 75 | _ε_11=15.21 | _ε_11=14.32 | _ε_11=11.7 |
| **SiC** | 182 | 7.10 | 6.01 | 6.552 |
| **AlP** | 1327 | 10.33 | 6.94 | 7.54 |
| **BN** | 17 | _ε_11=4.75 | _ε_11=3.72 | _ε_11=5.06 |
| **BP** | 1312 | 9.03 | 7.94 | 11.0 |
| **GaP** | 1393 | 13.22 | 8.33 | 11.11 |
| **AlSb** | 1408 | 12.27 | 9.87 | 12.04 |
| **ZnS** | 1702 | 9.39 | 4.8 | 8.0 |
| **CdTe** | 23 | 19.59 | 6.54 | 10.6 |
| **HgTe** | 8041 | _ε_11=29.44 | _ε_11=11.22 | _ε_11=20 |
| **ZnSiP2** | 2376 | _ε_11=12.44 | _ε_11=8.56 | _ε_11=11.15 |
| **ZnGeP2** | 2355 | _ε_11=14.75 | _ε_11=9.02 | _ε_11=15 |
| **MAE** | - | 2.46 | 2.78 | - |

![](RackMultipart20200915-4-htfc92_html_cbd168cc1e1dc7fa.png)

![](RackMultipart20200915-4-htfc92_html_128a5621eb654f54.png)

_Fig. Dielectric tensor (3x3 shape) distributions of the dataset considered._

**1.3.10 Piezoelectric tensor**

The piezoelectric effect is a reversible process where materials exhibit electrical polarization resulting from an applied mechanical stress, or conversely, a strain due to an applied electric field. Common applications for piezoelectricity include medical applications, energy harvesting devices, actuators, sensors, motors, atomic force microscopes, and high voltage power generation. PZ responses can be measured under constant strain, giving the piezoelectric stress tensor or constant stress, giving the piezoelectric strain tensor . The PZ data is obtained using DFPT method discussed above. Piezoelectric tensor is a 3x3 matrix. A distribution of the PZ tensor components is shown below:

![](RackMultipart20200915-4-htfc92_html_5d3e3d0d39021471.png)

![](RackMultipart20200915-4-htfc92_html_7eda49dfb35924a9.png)

_Fig. Piezoelectric tensor (3x6 shape) distributions of the dataset considered._

_Table Comparison of piezoelectric coefficient max(e __ij__ ) data for experiments and DFT. We take average values for the cases where the experimental data are in a range._

| **Mats.** | **JID** | **Max(eij)** | **DFT** | **Reference** |
| --- | --- | --- | --- | --- |
| **BN** | 57695 | 1.55 | 1.15 | 9 |
| **AlN** | 39 | 1.39-1.55 | 1.39 | 10-12 |
| **ZnS** | 7648 | 0.38 | 0.13 | 9 |
| **ZnSe** | 8047 | 0.35 | 0.06 | 9 |
| **SiO2** | 41 | 0.171 | 0.16 | 13,14 |
| **BaTiO3** | 110 | 1.94-5.48 | 4.13 | 15,1617,18 |
| **LiNbO3** | 1240 | 1.33 | 1.59 | 19 |
| **GaSb** | 35711 | -0.07 | -0.102 | 9 |
| **PbTiO3** | 3450 | 3.35-5.0 | 3.96 | 20 |
| **GaN** | 30 | 0.73 | 0.47 | 21 |
| **InN** | 1180 | 0.97 | 0.90 | 22 |
| **AlP** | 1327 | -0.1 | 0.004 | 9 |
| **AlAs** | 1372 | -0.16 | 0 | 9 |
| **AlSb** | 1408 | -0.13 | 0.06 | 9 |
| **ZnO** | 1195 | 1.00-1.55,0.89 | 1.10 | 23 |
| **BeO** | 20778 | 0.1 | 0.22 | 23 |
| **MAD** |
 |
 | 0.21 |
 |

**1.3.11 Infrared intensity**

The infrared intensity is important for thermal-imaging, infrared-astronomy, food-quality control. Infrared frequencies are classified in three categories: far (30-400 cm-1), mid (400-4000 cm-1) and near (4000-14000 cm-1) IR frequencies. The IR intensity is calculated obtained from the gamma-point phonon data used in the DFPT calculations.

_Table Comparison of experimental and DFPT IR frequencies (cm __-1__ )._

| **Mats.** | **JID** | **DFPT** | **Experiment** |
| --- | --- | --- | --- |
| **ZnO** | 1195 | 379, 410 | 389,4131 |
| **AlN** | 39 | 600, 653 | 620,6692 |
| **GaN** | 30 | 532 | 5313,4 |
| **SnS** | 1109 | 93, 144, 178, 214 | 99, 145, 178, 220 5 |
| **SnSe** | 299 | 72.6, 98.44, 125.01, 160.3 | 80,96,123, 150 5 |
| **KMgF3** | 20882 | 160.0, 287.0, 470.8 | 168, 299, 4586 |
| **LiNbO3** | 1240 | 145.0, 216.6 | 160, 2207 |
| **GeS** | 2169 | 106.4, 196.2, 236.5, 253.0, 276.9 | 118, 201,238,258, 2808 |
| **MAD** |
 |
 | 8.36 |

![](RackMultipart20200915-4-htfc92_html_33f80effe0a17282.png)

_Fig. Analysis of the IR-data. a) IR peaks for all the materials in the database, b) Comparison of finite-difference (FDM) and DFPT phonon frequencies for conventional cells, c) space-group distributions of materials with at least one far (blue) and mid (green) IR peaks, d) Venn-diagram of the chemistry of materials containing chalcogenides, halides. or pnictides. materials, e) Minimum frequency vs density of the system, f) dimensionality analysis of the far-IR materials_.

**1.3.12 Elastic tensor**

We use finite difference method on conventional cells of systems for ET calculations. For bulk material, the compliance tensor can be obtained by:

Now, several other elastic properties calculated from Cij and sij. Some of the important properties are given below:

KV = ((C11+C22+C33) + 2(C12+C23+C31))/9

GV = ((C11+C22+C33) − (C12+C23 + C31) +3 (C44+C55+C66))/15

KR = ((s11+s22+s33) + 2(s12+s23+s31))-1

GR = 15(4(s11+s22+s33) - 4(s12+s23+s31) + 3(s44+s55+s66))-1

KVRH =(KV+KR)/2

GVRH =(GV+GR)/2

ν = (3KVRH − 2GVRH)/(6KVRH+2GVRH))

Here KV and GV are Voigt bulk and shear modulus, and KR and GR Reuss-bulk and shear modulus respectively. The homogenous Poisson ratio is calculated as ν. The EC data can be also used to predict the ductile and brittle nature of materials with Pugh (Gv/Kv) and Pettifor criteria (C12-C44) . Materials with Pugh ratio value \&gt;0.571 and Pettifor criteria \&lt;0 should be brittle, while materials with Pugh ratio value \&lt;0.571 and Pettifor criteria \&gt;0 should be ductile.

For monolayer material calculations, the elastic tensor obtained from DFT code such as VASP, assumes periodic-boundary-condition (PBC). Therefore, cell vectors are used to calculate the area which again is used in computing stress. When dealing with the monolayer, an arbitrary vacuum padding is added in one of the direction (say z-direction). When computing EC we need to correct the output by eliminating the arbitrariness of the vacuum padding. We do that as a post-processing step by multiplying the _C__ij_ components by the length of the vacuum padding. Therefore, the units of EC turn into Nm-1 from Nm-2. For example, in order to calculate C11 (stress computed in x direction), the area is computed using normal of y and z-vectors. Obviously, the z-vector is arbitrary, so if we multiply the output by z-vector magnitude we get rid of the arbitrariness of z and also get C11 in Nm-1. As shown in Fig. 1, the z-vector magnitude is the z-simulation. The above discussion can also be expressed as the following:

![](RackMultipart20200915-4-htfc92_html_b662f2a41fb4d4e7.png)

![](RackMultipart20200915-4-htfc92_html_2939c22591a2a74c.png)

_Fig. Elastic constant distribution for 3D (magenta), 2D (green), 1D (blue) and 0D (red) materials._

_Table. Comparison of bulk modulus, K__V_ _(GPa), from vdW-DF-optB88 (OPT) and experiments. The experimental data are however not data corrected for zero-point energy effects, which would lead to a slight increase__of the values._

| **Material** | **JVASP#** | **OPT** | **Expt.** |
 | **Material** | **JVASP#** | **OPT** | **Expt.** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Cu** | 14648 | 141.4 | 142 |
 | **V** | 1041 | 183.4 | 161.9 |
| **C (diamond)** | 91 | 437.4 | 443 |
 | **Fe** | 882 | 193 | 168.3 |
| **Si** | 1002 | 87.3 | 99.2 |
 | **Ni** | 14630 | 200.4 | 186 |
| **Ge** | 890 | 58.1 | 75.8 |
 | **Nb** | 934 | 176 | 170.2 |
| **Ag** | 813 | 100.3 | 109 |
 | **Mo** | 925 | 262 | 272.5 |
| **Pd** | 14644 | 176 | 195 |
 | **Ta** | 14750 | 199 | 200 |
| **Rh** | 14817 | 260.8 | 269 |
 | **W** | 14830 | 305.2 | 323.2 |
| **Li** | 913 | 13.9 | 13.3 |
 | **Ir** | 901 | 348 | 355 |
| **Na** | 25140 | 7.7 | 7.5 |
 | **Pt** | 972 | 251.6 | 278.3 |
| **K** | 14800 | 3.9 | 3.7 |
 | **Au** | 825 | 148 | 173.2 |
| **Rb** | 978 | 3.1 | 2.9 |
 | **Pb** | 961 | 42.6 | 46.8 |
| **Ca** | 846 | 17.7 | 18.4 |
 | **LiCl** | 23864 | 35.5 | 35.4 |
| **Sr** | 21208 | 12.5 | 12.4 |
 | **NaCl** | 23862 | 27.7 | 26.6 |
| **Ba** | 831 | 9.9 | 9.3 |
 | **NaF** | 20326 | 53.7 | 51.4 |
| **Al** | 816 | 70 | 79.4 |
 | **MgO** | 116 | 160.7 | 165 |
| **LiF** | 1130 | 73.9 | 69.8 |
 | **SiC** | 182 | 213.3 | 225 |
| **TiO2-anatase** | 314 | 196 | 191.9 |
 | **GaAs** | 1174 | 62 | 75.6 |
| **TiO2-rutile** | 10036 | 226.3 | 243.5 |
 | **P (black)** | 7818 | 41 | 36 |
| MAE (GPa): | 8.51 |
 |
 |
 |
 |
 |
 |
 |

![](RackMultipart20200915-4-htfc92_html_ed6769ebaaaadeeb.png)

_Fig. Periodic table trend for high bulk modulus material constituents. The bulk moduli of all the materials were projected on individual elements and their average contribution is shown. The colorbar is in the unit of GPa. A similar trend was found for shear modulus._

From DFPT and finite-difference calculations, force constants (mainly at gamma point in our case) can be calculated. Using this information, we can predict the phonon density of states and bandstructure assuming the conventional cell used in the calculation is large enough (at least 11 Angstrom). Note that for smaller lengths the phonon bandstructure and density of states can show unreasonable frequencies at off and away from gamma points.

![](RackMultipart20200915-4-htfc92_html_525beddac7fb875a.png)

![](RackMultipart20200915-4-htfc92_html_abf574c4094a878e.png)

**1.3.13 Thermoelectric properties:**

Thermoelectrics[1-4](#_ENREF_1) are materials that can convert a temperature gradient into electric voltage, or vice-versa. Themoelectrics can be used to regenerate electricity from waste heat[5](#_ENREF_5), refrigeration[6](#_ENREF_6) and several other space-technology applications7,8. The search for efficient thermoelectric materials is an area of intense research due the potential of converting waste heat into electrical power, and therefore improving energy efficiency and reducing fossil fuel usage. The figure of merit for thermoelectric materials is the dimensionless quantity _zT_:

where _S_, , _k __e__ , k__l_, and _T_ are the Seebeck coefficient, electrical conductivity, electronic part of thermal conductivity, lattice thermal conductivity, and temperature, respectively. The numerator, , is referred to as the power-factor. To achieve a high _zT_, a material should have a high-power factor and low thermal conductivity.

_Table : Chemical formula, experimental Seebeck value (_μV/K)_, DFT value, JARVIS-ID, doping concentration, doping type, temperature, space-group and reference data for the experimental vs DFT comparisons._

| **Formula** | **Exp** | **DFT** | **JID** | **Dop.conc.**** (/cm3) **|** type **|** T****(K)** | **Spg.** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **Bi2Te3** | 116 | 124.7357 | JVASP-25 | 7.78E+19 | p | 420 | 166 |
| **Bi2Se3** | -70 | -136.9 | JVASP-1067 | -2.20E+19 | n | 420 | 166 |
| **CuInTe2** | 254 | 203.5364 | JVASP-3495 | 1.60E+19 | P | 300 | 122 |
| **CuGaTe2** | 380 | 448.9839 | JVASP-2295 | 1.00E+18 | p | 300 | 122 |
| **AgTlTe** | 550 | 721 | JVASP-9744 | 1.00E+17 | p | 320 | 62 |
| **ErNiSb** | 258 | 268.71 | JVASP-1903 | 1.42E+19 | p | 335 | 216 |
| **Cu2ZnSnSe4** | -26.02 | -23.98 | JVASP-17430 | 1.00E+18 | p | 293 | 121 |
| **CoNbSn** | -69 | -2.22 | JVASP-18668 | 5.97E+16 | p | 318 | 216 |
| **AlFe2V** | -107 | -32.3911 | JVASP-15637 | 5.00E+20 | n | 300 | 225 |
| **CoSbZr** | -62 | -43.9953 | JVASP-18207 | 2.72E+20 | n | 335 | 216 |
| **SnSe** | 586 | 674.7 | JVASP-299 | 3.16E+17 | p | 523 | 62 |
| **SnTe** | 103 | 111.855 | JVASP-7860 | 1.00E+21 | p | 817 | 225 |
| **Cu2Se** | 258 | 148.2624 | JVASP-18192 | 2.00E+21 | p | 900 | 216 |
| **Mg2Sn** | -71.5 | -91.3387 | JVASP-14507 | -2.00E+19 | n | 400 | 225 |

![](RackMultipart20200915-4-htfc92_html_55493afc10ad0268.png)

_Fig. Flow chart associated with the data-driven thermoelectric materials design._

![](RackMultipart20200915-4-htfc92_html_8272ddda79482333.png)

![](RackMultipart20200915-4-htfc92_html_3c039f3c9eb6880c.png)

_Fig. A brief overview of the thermoelectric data for periodic bulk materials. Figure a) -f) have been computed on all 3D materials, while g-i) display properties only computed on the set of bulk screened materials (bandgap\&gt;0.05 eV, n-type and p-type PF\&gt;1000 (μW/(mK) __2__ ) at 600 K and 10 __20_ _cm__ -3_ _doping concentration). a) n-type Seebeck coefficient distribution, b) n-type power factor of materials, c) n-type electrical conducitivity plotted against the absolute values of Seebeck-coefficient with color-coded power-factor and size of the dots proportional to bandgaps, d) p-type Seebeck coefficient distribution, e) p-type power factor of materials, f) p-type power factor plotted against the bandgaps, g) predicted dimensionality distribution of screened materials, h) space-group distribution of the screened materials, i) power-factor vs Seebeck distribution of the screened materials._

![](RackMultipart20200915-4-htfc92_html_e934d6f35ac333a8.png)

_Fig. Periodic table trends of high-power factor materials. The elements in a material were assigned 1 or 0 if the material has high or low n-type power-factor (\&gt;1000 (μW/(mK) __2__ )). Then the probability of finding the element in a high power-factor material was calculated._

**1.3.14 Wannier tight binding Hamiltonians**

Wannier functions (WF) were first introduced in 1937, and have proven to be a powerful tool in the investigation of solid-state phenomenon such as polarization, topology, and magnetization. WTBH is not necessarily a material properties but can be useful in calculating several material properties. Mathematically, WFs are a complete orthonormalized basis set that act as a bridge between a delocalized plane wave representation commonly used in electronic structure calculations and a localized atomic orbital basis that more naturally describes chemical bonds. One of the most common ways of obtaining Wannier tight-bonding Hamiltonians (WTBH) is by using the Wannier90 software package to generate maximally localized Wannier functions, based on underlying density functional theory (DFT) calculations.

![](RackMultipart20200915-4-htfc92_html_43ca0b015464fec3.png)

_Fig. Workflow showing the Wannierization from using the DFT calculations._

Wannier bands compared with explicit DFT bands to assess the accuracy of the Wannier tight-binding Hamiltonians. For example:

![](RackMultipart20200915-4-htfc92_html_54da7f7e26ae7b75.png)

_Fig. Comparison of DFT and WTB bandstructures for Bi __2__ Se __3__. a-b) on dense k-grid, c-d) high-symmetry Brillouin zone points._

![](RackMultipart20200915-4-htfc92_html_23ed3f4a12be1610.png)

_Fig. DFT-TB maximum difference (μ) distribution for all the Wannier Tight-binding Hamiltonians (WTBHs). A) on a regular k-point grid, b) on high-symmetry k-points._

**1.3.15 Scanning tunneling microscopy images**

Since the invention of the scanning tunneling microscope (STM), this technique has become an essential tool for characterizing material surfaces and adsorbates. In addition to providing atomic insights, STM has been proven useful for characterizing the electronic structure, shapes of molecular orbitals, and vibrational and magnetic excitations. It can also be used for manipulating adsorbates and adatoms, and for catalysis and quantum information processing applications. Few examples of computational STM is given below:

![](RackMultipart20200915-4-htfc92_html_1d7f497b58ddbc7a.png)

_Fig. STM images for which we qualitatively compared computational positive bias images._ We compare computational STM images with that of experiments for graphene16, 2H-MoS217, 2H-NbSe218, 2H-WSe219, 1T&#39;-WTe220, FeSe21, black-P22,23, SnSe24, Bismuth surface57,58. _The JVASP identifiers can be used to visualize the detailed webpage for each material._

![](RackMultipart20200915-4-htfc92_html_98f8a44ea3ed2f24.png)

**1.3.16 Electric field gradients**

The key parameters used to define Nuclear Quadrupole Resonance (NQR) spectral lines are the quadrupole coupling constant Q = and the asymmetry parameter

,

where _e_ is electric charge, _h_ is Planck&#39;s constant, and _Q_ is the nuclear quadrupole moment; _V __ii_are the principal components of the diagonalized EFG tensor, defined as the second derivative in Cartesian coordinates of the Coulomb potential at the nucleus position. By construction, the EFG, _V__ ii_ is a traceless tensor. The coordinate system, in accordance with the convention used by experimentalists is chosen so that , which forces . Note that if the point group of the site in question is cubic, then by symmetry all components are zero; if it is tetragonal or hexagonal, then __= 0, but _V_zz 0.

Nuclear resonance spectroscopies, such as Nuclear Magnetic Resonance (NMR) and Nuclear Quadrupole Resonance (NQR) are extremely valuable as sensitive probes of the local electronic structure in solids. They are considered to be the golden standard in addressing such disparate physical properties of materials such as magnetism, charge ordering, structural distortion, valence skipping, superconductivity, and many others. Shifts of the NMR/NQR spectral lines yield information of uniform susceptibilities, while their relaxation time informs about the local susceptibility at a nucleus. The very low excitation energies of nuclear resonances (even on the scale of the superconducting gap) provide a window into dynamical effects. Further, nuclear resonances have a plethora of practical applications, ranging from medical (MRI) and pharmaceutical4-6 to detecting prohibited substances and combatting terrorism.

The difference between NMR and NQR is that in NMR, the separation of the nuclear levels is predominantly affected through external or internal magnetic fields, while in NQR it comes from the interaction of the nuclear quadrupolar moment with the gradients of the static electric field at the nucleus (Electric Field Gradients, EFG). EFG can also be obtained from Mossbauer spectroscopy, however, NMR and NQR have recently become more commonly used especially for crystalline solids.

_Table Comparison of current density functional (J-DFT) predictions with experimental (Exp) and previously (Prev.-DFT) reported Electric Field Gradient, V__zz_ _(10__21 __Vm__ -2__) data. The MAD (Mean Absolute Deviation), and MAPD (Mean Absolute Percentage Difference) values are calculated for the whole data. Details of each material are provided at its corresponding webpage. Please look into the references (and references therein) for experimental and previously calculated data._

| **Material** | **JID** | **Atom** | **|V**** zz (Exp) ****|** | **|V**** zz (J-DFT) ****|**
 | **|V**** zz (Prev.-DFT) ****|**
 | **|Δ|** | **|Δ%|** |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Cl2 | 855 | Cl | 55.1858 | 52.85 | 54.2343 | 2.33 | 4.22 |
| Br2 | 840 | Br | 95.6958 | 88.86 | 94.4443 | 6.83 | 7.14 |
| I2 | 895 | I | 113.0058 | 108.70 | 119.0143 | 4.30 | 3.81 |
| Be | 25056 | Be | 0.04459 | 0.072 | 0.0660 | 0.028 | 63.64 |
| Mg | 14840 | Mg | 0.04859 | 0.079 | 0.0460 | 0.031 | 64.58 |
| Sc | 996 | Sc | 0.3859 | 1.78 | 0.9660 | 1.40 | 368.4 |
| Ti | 14815 | Ti | 1.6159,61 | 1.64 | 1.7560 | 0.03 | 1.86 |
| Co | 858 | Co | 2.959 | 0.52 | 0.2960 | 2.38 | 82.06 |
| Zn | 1056 | Zn | 3.4859 | 5.62 | 4.2960 | 2.14 | 61.50 |
| Zr | 14612 | Zr | 4.4059 | 3.50 | 4.1460 | 0.90 | 20.45 |
| Tc | 1020 | Tc | 1.8359 | 1.67 | 1.7460 | 0.16 | 8.74 |
| Ru | 987 | Ru | 0.9759 | 1.52 | 1.6260 | 0.55 | 56.70 |
| Cd | 14832 | Cd | 6.5059 | 7.56 | 8.1360 | 1.06 | 16.31 |
| La | 910 | La | 1.6259 | 2.24 | 0.9160 | 0.62 | 38.27 |
| Hf | 14590 | Hf | 7.3359 | 8.87 | 8.1260 | 1.54 | 21.01 |
| Re | 981 | Re | 5.1259 | 6.14 | 6.4960 | 1.02 | 19.92 |
| Os | 952 | Os | 4.1659 | 6.00 | 7.0260 | 1.84 | 44.23 |
| BI3 | 3630 | I | 71.2962 | 68.98 | - | 2.31 | 3.24 |
| CF3I | 32512 | I | 124.3463 | 123.22 | - | 1.12 | 0.90 |
| CIN | 5758 | I | 157.2164 | 151.0 | - | 6.21 | 3.95 |
| NaNO2 | 1429 | Na | 0.43865 | 0.552 | 0.57566 | 0.114 | 26.03 |
| NaNO2 | 1429 | N | 11.1065 | 12.194 | 11.77266 | 1.094 | 9.86 |
| Cu2O | 1216 | Cu | 9.8067 | 6.47 | 6.76566 | 3.33 | 33.98 |
| TiO2 | 10036 | Ti | 2.2168 | 2.098 | 2.26966 | 0.112 | 5.07 |
| TiO2 | 10036 | O | 2.3868 | 2.21 | 2.23566 | 0.17 | 7.14 |
| SrTiO3 | 8082 | O | 1.6269 | 1.24 | 1.0070 | 0.38 | 23.46 |
| BaTiO3 | 8029 | O | 2.4669 | 3.56 | 2.3570 | 1.10 | 44.72 |
| Li3N | 1375 | N | 1.0471 | 1.25 | 1.0966 | 0.21 | 20.19 |
| Li3N | 1375 | Li(2c) | 0.3071 | 0.225 | 0.29166 | 0.075 | 25.00 |
| Li3N | 1375 | Li(1b) | 0.6071 | 0.50 | 0.61666 | 0.144 | 24.00 |
| FeSi | 8178 | Fe | 4.4572,73 | 4.84 | 4.9243 | 0.39 | 8.76 |
| FeS2(marcasite) | 2142 | Fe | 3.074 | 2.93 | 3.2143 | 0.07 | 2.33 |
| FeS2(pyrite) | 9117 | Fe | 3.6674 | 3.51 | 3.4043 | 0.15 | 4.10 |
| 2H-MoS2 | 54 | Mo | 7.0939,40 | 7.70 | - | 0.61 | 8.60 |
| 2H-MoS2 | 54 | S | 5.5439,40 | 5.33 | - | 0.21 | 3.80 |
| 2H-WS2 | 72 | S | 4.8239,40 | 4.53 | - | 0.29 | 6.22 |
| CaGa2 | 16464 | Ga | 4.4475 | 3.55 | 3.7776 | 0.89 | 20.05 |
| SrGa2 | 14853 | Ga | 5.2275 | 2.54 | 4.1376 | 2.68 | 51.34 |
| BaGa2 | 19628 | Ga | 4.4875 | 5.10 | 4.3876 | 0.62 | 13.84 |
| NaGa4 | 14728 | Ga(e) | 6.4976 | 5.20 | 6.1876 | 1.29 | 19.88 |
| NaGa4 | 14728 | Ga(d) | 4.6476 | 4.33 | 4.4476 | 0.31 | 6.68 |
| CaGa4 | 20533 | Ga(e) | 2.8976 | 2.67 | 2.8076 | 0.22 | 7.61 |
| CaGa4 | 20533 | Ga(d) | 4.8776 | 4.99 | 4.7376 | 0.12 | 2.46 |
| SrGa4 | 20206 | Ga(e) | 2.5176 | 1.67 | 2.2476 | 0.84 | 33.47 |
| SrGa4 | 20206 | Ga(d) | 5.9576 | 5.31 | 5.6476 | 0.64 | 10.76 |
| TaP | 79643 | Ta | 3.0077 | 2.5 | 3.5477 | 0.50 | 16.67 |
| UAs2 | 19797 | U | 15.078 | 9.7 | 13.0357 | 5.3 | 35.3 |
| **MAD** |
 |
 |
 |
 |
 | 1.17 | - |
| **MAPD** |
 |
 |
 |
 |
 | - | 28.90 |

![](RackMultipart20200915-4-htfc92_html_a6140f133a7ab9bc.png)

_Fig. Histogram (above) and intensity map (below) of the probability are shown as a function of VZZ and . The logarithmic z-scale allows for the strong preponderance of small VZZ and small  values to be clearly seen, while at the same time revealing the bell-shaped distribution in VZZ for all  values. Values where Vzz=0, and therefore  is ill-defined, have been excluded._

**1.3.17 DFT convergence parameters**

Although density functional theory is exact in theory, its implementation requires several approximations such as the choice of basis-set, exchange-correlation functional, mesh-size for Brillouin zone (BZ) integration and plane-wave cut-off for plane-wave basis. These parameters need to be converged prior to geometric optimization and property predictions. Such convergences are performed for k-point and plane wave cut-off only in JARVIS-DFT leading to high quality of the data.

![](RackMultipart20200915-4-htfc92_html_38dafc26ed259ad3.png)

_Fig. An example of convergence in a) cut-off and b) k-points for Bi __2__ Se__3_

![](RackMultipart20200915-4-htfc92_html_b3eea6373a7f3f79.jpg)

_Fig. Histogram for ranges of the number of k-points. The &#39;energy per cell&#39; and &#39;energy per atom&#39; methods are denoted by &#39;EPC&#39; and &#39;EPA&#39; respectively. EPA values are derived from EPC. Fig. a) shows the length-based k-points distribution of all the materials in the database using EPC, b) the per reciprocal-atom-based distribution of k-points using EPC, c) length-based k-points distribution of all the materials in the database using EPA, d) the per reciprocal atom-based distribution of k-points using EPA. Similar distributions for non-metals are shown in e-h and for metals in i-l. Metals are in general observed to require more k-points._

**1.4 Data Quality assessment table**

_Table. Mean Absolute Error (MAE) for JARVIS-DFT data with respect to available experimental data for various material properties._

| **Property** | **#Materials** | **MAE** | **Typical range** |
| --- | --- | --- | --- |
| **Formation energy (eV/atom)** | 1317 | 0.128 | -4 to 2 |
| **OptB88vdW-bandgaps (eV)** | 54 | 1.33 | 0 to 10 |
| **TBmBJ-bandgaps (eV)** | 54 | 0.51 | 0 to 10 |
| **Bulk modulus (GPa)** | 21 | 5.75 | 0 to 250 |
| **Electronic (𝜀11) OPT** | 28 | 3.2 | 0 to 60 |
| **Electronic (𝜀11) MBJ** | 28 | 2.62 | 0 to 60 |
| **Solar-eff. (SLME) (%) (MBJ)** | 5 | 6.55 | 0 to 33.7 |
| **Max. piezoelectric strain coeff (Cm-2)** | 16 | 0.21 | 0 to 2 |
| **Dielectric constant (𝜀11) (DFPT)** | 16 | 2.46 | 0 to 60 |
| **Seebeck coefficient (μV/K)** | 14 | 54.7 | -600 to 600 |
| **Electric field gradient Vzz (1021Vm-2)** | 37 | 1.17 | 0 to 100 |
| **IR mode (cm-1)** | 8 | 8.36 | 0 to 4000 |

**1.5 Statistics**

_Table. A brief summary of datasets available in the JARVIS-DFT._

| **Material classes** | **Numbers** |
| --- | --- |
| **3D-bulk** | 33482 |
| **2D-bulk** | 2293 |
| **1D-bulk** | 235 |
| **0D-bulk** | 413 |
| **2D-monolayer** | 1105 |
| **2D-bilayer** | 102 |
| **Molecules** | 12 |
| **Heterostructure** | 3 |
| **Total DFT calculated systems** | 37646 |

_Table. A brief summary of functionals used in optimizing crystal geometry in the JARVIS-DFT._

| **Functionals** | **Numbers** |
| --- | --- |
| **vdW-DF-OptB88 (OPT)** | 37646 |
| **vdW-DF-OptB86b (MK)** | 109 |
| **vdW-DF-OptPBE (OR)** | 111 |
| **PBE** | 99 |
| **LDA** | 92 |

_Table. A brief summary of material-properties available in the JARVIS-DFT. The database is continuously expanding_.

| **JARVIS-DFT DATA Property** | **Numbers** |
| --- | --- |
| **Optimized crystal-structure (OPT)** | 37646 |
| **Formation-energy (OPT)** | 37646 |
| **Bandgap (OPT)** | 37646 |
| **Exfoliation energy (OPT)** | 819 |
| **Bandgap (TBmBJ)** | 15655 |
| **Bandgap (HSE06)** | 40 |
| **Bandgap (PBE0)** | 40 |
| **Bandgap (G****0 ****W**** 0****)** | 15 |
| **Bandgap (DMFT)** | 11 |
| **Frequency dependent dielectric tensor (OPT)** | 34045 |
| **Frequency dependent dielectric tensor (TBmBJ)** | 15655 |
| **Elastic-constants (OPT)** | 15500 |
| **Finite-difference phonons at Г-point (OPT)** | 15500 |
| **Work-function, electron-affinity (OPT)** | 1105 |
| **Theoretical solar-cell efficiency (SLME) (TBmBJ)** | 5097 |
| **Topological spin-orbit spillage (PBE+SOC)** | 11500 |
| **Wannier tight-binding Hamiltonians (PBE+SOC)** | 1771 |
| **Seebeck coefficient (OPT, BoltzTrap)** | 22190 |
| **Power factor (OPT, BoltzTrap)** | 22190 |
| **Effective mass (OPT, BoltzTrap)** | 22190 |
| **Magnetic moment (OPT)** | 37528 |
| **Piezoelectric constant (OPT, DFPT)** | 5015 |
| **Dielectric tensor (OPT, DFPT)** | 5015 |
| **Infrared intensity (OPT, DFPT)** | 5015 |
| **DFPT phonons at Г-point (OPT)** | 5015 |
| **Electric field gradient (OPT)** | 15187 |
| **Non-resonant Raman intensity (OPT, DFPT)** | 250 |
| **Scanning tunneling microscopy images (PBE+SOC)** | 770 |

**1.6 References**

1. _JARVIS: An Integrated Infrastructure for Data-driven Materials Design, arXiv:2007.01831 (2020)._
2. _High-throughput Identification and Characterization of Two-dimensional Materials using Density functional theory, Scientific Reports 7, 5179 (2017)._
3. _Computational Screening of High-performance Optoelectronic Materials using OptB88vdW and TBmBJ Formalisms, Scientific Data 5, 180082 (2018)_
4. _Elastic properties of bulk and low-dimensional materials using Van der Waals density functional, Phys. Rev. B, 98, 014107 (2018)._
5. _High-throughput Discovery of Topologically Non-trivial Materials using Spin-orbit Spillage, Scientific Reports 9, 8534 (2019)_
6. _Convergence and machine learning predictions of Monkhorst-Pack k-points and plane-wave cut-off in high-throughput DFT calculations, Comp. Mat. Sci. 161, 300 (2019)._
7. _High-throughput Density Functional Perturbation Theory and Machine Learning Predictions of Infrared, Piezoelectric and Dielectric Responses,npj Computational Materials, 6, 64 (2020)._
8. _Computational Search for Magnetic and Non-magnetic 2D Topological Materials using Unified Spin-orbit Spillage Screening, npj Computational Materials, 6, 49 (2020)._
9. _Accelerated Discovery of Efficient Solar-cell Materials using Quantum and Machine-learning Methods, Chem. Mater., 31, 15, 5900 (2019)._
10. _Data-driven Discovery of 3D and 2D Thermoelectric Materials, J. Phys.: Condens. Matter 32, 475501 (2020)._
11. _Efficient Computational Design of 2D van der Waals Heterostructures: Band-Alignment, Lattice-Mismatch, Web-app Generation and Machine-learning, submitted._
12. _Density Functional Theory based Electric Field Gradient Database, arXiv:2005.09255._
13. _Database of Wannier Tight-binding Hamiltonians using High-throughput Density Functional Theory, arXiv:2007.01205._
14. _Density Functional Theory and Deep-learning to Accelerate Data Analytics in Scanning Tunneling Microscopy, arXiv:1912.09027._

**2 JARVIS-FF**

**2.1 Overview**

JARVIS-FF is a repository of classical force-field/potential calculation data intended to help users select the most appropriate force-field for a specific application. Many classical force-fields are developed for a particular set of properties (such as energies), and may not have been tested for properties not included in training (such as elastic constants, or defect formation energies). JARVIS-FF provides an automatic framework to consistently calculate and compare basic properties, such as the bulk modulus, defect formation energies, phonons, etc. that may be critical for specific molecular-dynamics simulations. JARVIS-FF relies on DFT and experimental data to evaluate accuracy.

![](RackMultipart20200915-4-htfc92_html_c8e013e969479ff2.png)

  1. **Methodology**

We started by downloading all the available potentials from the NIST interatomic potential repository (IPR) and from LAMMPS itself (15May2015 version). For each element having at least a potential, we downloaded all the corresponding crystal structures from the Materials Project (MP) or JARVIS-DFT database. We also downloaded all the energetics and mechanical properties data from and stored them in a separate database for a later comparison with the classical results (LAMMPS calculations at T=0 K). The high-throughput setting of LAMMPS jobs was done using JARVIS-tools. In our runs we used 10-06 as strain, 10-10 eV/Å for force convergence during the minimization to optimize the structure and 1000 maximum iteration for structure optimization. These are generalized computational set-up parameters, and the energetics and elastic constant data may or may not depend on them. We tested strain parameters for a range of values (10-04, 10-06 and 10-08) but obviously evaluating such set of parameters for all the calculations was too extensive a work and was not carried out here. The relaxed structure was also stored along with the above files for later use such as for performing defect, phonon or other similar calculations.

After the minimization, the crystal structure is stored in LAMMPS data-format and JSON format. Using this JSON file, unique Wyckoff positions were identified and deleted to represent vacancy-structures. The multiplicity of the Wyckoff positions is also recorded. After the defect structure generation, the LAMMPS energy minimization is carried out. In a subsequent run, we calculate the chemical potential of the defect element using the specific force-field.

The data for the vacancy structure, chemical potential of element and perfect structure energy were used to calculate the defect formation energies. The most stable prototype for chemical potential calculation was determined using the energy above convex hull data from DFT. The defect structures were required to be at least 1.5 nm long in the x, y and z directions to avoid self-interactions with the periodic images of the simulation cell. Similar to the defect structures, distinct surfaces were created up to 3 Miller indices with the relaxed structure stored in the JSON file. A generic code for generating defect and surface structures is given at our github page. We enforce the surfaces to be at least 2.5 nm thick and with 2.5 nm vacuum in the simulation box. The 2.5 nm vacuum is used to ensure no self-interaction and the thickness is used to mimic actual experimental surface. Using the energies of perfect bulk and surface structures, surface energies for a specific plane are calculated. We should point out that only unreconstructed surfaces without any surface-segregation effects are computed, as our high-throughput approach does not allow for taking into account specific, element dependent reconstructions yet.

Phonons were obtained by making an interface of JARVIS-FF with Phonopy package at 0 K. For deformed-structures, constant volume ensemble was used. The deformed structures were taken of at least 1.5 nm size in all directions.

![](RackMultipart20200915-4-htfc92_html_a87caa4e69f4e47.png)

_Fig. Flow-chart showing the processes involved in creating the JARVIS-FF database._

_Table. A summary of various types of force-fields available in the JARVIS-FF._

| **Force-fields** | **Numbers** |
| --- | --- |
| **EAM** | 92 |
| **Tersoff** | 9 |
| **ReaxFF** | 5 |
| **COMB** | 6 |
| **AIREBO** | 2 |
| **MEAM** | 1 |
| **EIM** | 1 |

  1. **Property details and assessment**

Using jarvis.core.Atoms class several atomistic properties such as lattice parameters, density, packing fraction etc. can be calculated. The optimized lattice parameters generally compare well with DFT data except for the FFs where a particular phase was not trained during FF-fitting. Similarly, energetics in terms of convex hull plot is compared between DFT and FF results. An example is shown below:

![](RackMultipart20200915-4-htfc92_html_eb2d372587730ea.jpg)

_Figure. Convex hull plot for Ni-Al system using a) DFT data, b) Force-field (Mishin Ni-Al potential). Filled blue circles show stable while the red squares show unstable structures. Here &quot;fu&quot; stands for formula unit._

Similarly, elastic tensor and derived properties were predicted using LAMMPS runs.

![](RackMultipart20200915-4-htfc92_html_d9889346cb9be914.png)

_Figure. Comparison of Voigt-bulk modulus data of force-fields and DFT data. The x-axis shows the bulk modulus for materials (for which FFs were available to us) from DFT calculations, while the y-axis shows the bulk modulus for corresponding materials using the force-field calculations._

Moreover, vacancy formation energies, surface energies can be compared between FF and DFT for FF quality assessment.

![](RackMultipart20200915-4-htfc92_html_208e8934e0f7fa61.png)

_Fig. Vacancy and surface energies (up to a maximum Miller index of 3) for single elements obtained with all the FF available for such materials are compared to available DFT data._

Phonon data from JARVIS-DFT and JARVIS-FF can be compared for a system to evaluate the phonon quality. However, it is important to note that in DFT there might be only conventional cell Gamma-point phonon data available whereas in JARVIS-FF we use supercell finite difference-based approach for obtaining phonon density of states and bandstructures. For a known stable material if the phonon bandstructure shows high negative values then it signifies the FF maynot be suitable to predict correct dynamical properties of the system.

  1. **References:**

_1. Evaluation and comparison of classical interatomic potentials through a user-friendly interactive web-interface, Scientific Data 4, 160125 (2017)._

_2. High-throughput assessment of vacancy formation and surface energies of materials using classical force-fields, J. Phys. Cond. Matt. 30, 395901(2018)._

**3 JARVIS-ML**

**3.1 Overview**

Machine learning has several applications in materials science and engineering, such as automating experimental data analysis, discovering new functional materials, optimizing known ones by accelerating conventional methods such as DFT, automating literature searches, discovering new physical equations, and efficient clustering of materials and their properties. There are several data types that can be used in ML such as scalar data (e.g., formation energies, bandgaps), vector/spectra data (e.g., density of states, dielectric function, charge density, X-ray diffraction patterns, etc.), image-based data (such as scanning tunneling microscopy and transmission electron microscopy images), and natural language processing-based data (such as scientific papers). In addition, ML can be applied on a variety of materials classes such as bulk crystals, molecules, proteins and free-surfaces.

**3.2 Methodology**

Currently, there are two types of data that are machine-learned in JARVIS-ML: discrete and image-based. The discrete target is obtained from the JARVIS-DFT database for 3D and 2D materials. There have been several descriptor developments as attempts to capture the complex chemical-structural information of a material. We compute CFID descriptors for most crystal structures in various databases. Many of these structures are non-unique but can still be used for pre-screening applications. The CFID can also be applied to other materials classes such as molecules, proteins, point defects, free surfaces, and heterostructures, which are currently ongoing projects. These descriptor datasets, along with JARVIS-DFT and other databases, act as input and outputs for machine learning algorithms. The CFID consists of 1557 descriptors for each material: 438 average chemical, 4 simulation-box-size, 378 radial charge-distribution, 100 radial distribution, 179 angle-distribution up to first neighbor, and another 179 for the second neighbor, 179 dihedral angle up to fist neighbor and 100 nearest neighbor descriptors.

![](RackMultipart20200915-4-htfc92_html_456fcd44ff109bc9.png)

_Fig. Figure showing different components of classical force-field inspired descriptors (CFID) for Si diamond structure. a) average radial-charge density distribution of constituent elements. b) total radial distribution function of the crystal structure, c) total angle distribution function up to the first-nearest neighbor, d) total dihedral-angle distribution up to the first-nearest neighbor, e) average chemical properties of constituent elements. The nearest neighbor distribution was obtained like the RDF._

**3.3 Quality assessment**

JARVIS-ML model accuracy is evaluated on the test-set (usually 10 %) representing newly computed and previously unseen DFT data for both regression and classifications models. Accuracy of regression and classification models are reported in terms of mean absolute error (MAE) and Receiver Operating Characteristic (ROC) Area Under Curve (AUC) metric respectively. A brief summary of regression and classification model accuracy results is given below in Table.

_Table. Performance of regression machine learning models in JARVIS-ML with JARVIS-DFT data using OptB88vdW (OPT) and TBmBJ (MBJ) with mean absolute error (MAE). The mean absolute deviation (MAD) of properties are also included._

| **Property** | **Training data** | **MAE** | **MAD** |
| --- | --- | --- | --- |
| **Formation energy (eV/atom)** | 24549 | 0.12 | 0.81 |
| **OPT bandgap (eV)** | 22404 | 0.32 | 1.05 |
| **MBJ bandgap (eV)** | 10499 | 0.44 | 1.60 |
| **Bulk mod., Kv (GPa)** | 10954 | 10.5 | 49.95 |
| **Shear mod., Gv (GPa)** | 10954 | 9.5 | 23.26 |
| **Refr. Index(x) (OPT)** | 12299 | 0.54 | 1.15 |
| **Refr. Index(x) (MBJ)** | 6628 | 0.45 | 1.03 |
| **IR mode (OPT) (cm-1)** | 3411 | 77.84 | 316.7 |
| **Max. Born eff. charge (OPT)(e)** | 3411 | 0.60 | 1.48 |
| **Plane-wave cutoff (OPT)(eV)** | 24549 | 85.0 | 370.6 |
| **K-point length (OPT)(Å)** | 24549 | 9.09 | 22.23 |
| **2D-Exfoliation energy(OPT) (eV/atom)** | 616 | 37.3 | 46.09 |

_Table. Performance of the classification machine learning models in JARVIS-ML with JARVIS-DFT data using OptB88vdW (OPT) and TBmBJ (MBJ) with Receiver Operating Characteristic (ROC) Area Under Curve (AUC) metric. Random guessing and perfect ROC AUC are 0.5 and 1 respectively._

| **Property** | **Number of datapoints** | **ROC AUC** |
| --- | --- | --- |
| **Metal/non-metal (OPT)** | 24549 | 0.95 |
| **Magnetic/Non-magnetic (OPT)** | 24549 | 0.96 |
| **High/low solar-cell efficiency (TBmBJ)** | 5097 | 0.90 |
| **High/low piezoelectric coeff** | 3411 | 0.86 |
| **High/low Dielectric** | 3411 | 0.93 |
| **High/low n-Seebeck coeff** | 21899 | 0.95 |
| **High/low n-power factor** | 21899 | 0.80 |
| **High/low p-Seebeck coeff** | 21899 | 0.96 |
| **High/low p-power factor** | 21899 | 0.82 |

![](RackMultipart20200915-4-htfc92_html_aebb479fb2b1710e.png)

_Fig. Performance and interpretation of formation energy ML model. a) learning curve, b) ML prediction on 10 % held samples, c) group wise feature importance of descriptors, d) comparison of model performance by incrementally adding various structural descriptors._

**3.4 References**

_1.     Machine learning with force-field inspired descriptors for materials: fast screening and mapping energy landscape, Phys. Rev. Mat., 2, 083801 (2018)._

_2.     Materials Science in the AI age: high-throughput library generation, machine learning and a pathway from correlations to the underpinning physics, MRS Comm., 9, 821 (2019)._

_3.     Enhancing materials property prediction by leveraging computational and experimental data using deep transfer learning, Nat. Comm., 10, 5316 (2019)._

_4.     Convergence and machine learning predictions of Monkhorst-Pack k-points and plane-wave cut-off in high-throughput DFT calculations, Comp. Mat. Sci. 161, 300 (2019)._

_5.   Accelerated Discovery of Efficient Solar-cell Materials using Quantum and Machine-learning Methods, Chem. Mater., 31, 15, 5900 (2019)._

_6.     Data-driven Discovery of 3D and 2D Thermoelectric Materials, arXiv:1906.06024._

_7.   High-throughput Density Functional Perturbation Theory and Machine Learning Predictions of Infrared, Piezoelectric and Dielectric Responses, npj Computational Materials, 6, 64 (2020)._

_8.  Density Functional Theory and Deep-learning to Accelerate Data Analytics in Scanning Tunneling Microscopy, arXiv:1912.09027._

65
