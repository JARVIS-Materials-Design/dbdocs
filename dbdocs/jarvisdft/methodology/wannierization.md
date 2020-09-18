We use Wannier90 to construct Maximally-Localized Wannier Functions (MLWF) based TB-Hamiltonians. The basic formalism of Wannierization is well-established. We briefly review some aspects here, interested readers can see longer discussions in 5,16. For a set of Bloch eigenvectors |├ ψ_(n,k) ⟩, a general set of WFs ├ |Rn⟩ (n=1...N) can be written as:
├ |Rn⟩=V/(2π)^3  ∫_BZ^BZ▒〖∑_(m=1)^N▒U_mn^((k)) |├ ψ_mk ⟩ e^(-ik.R) dk〗							
where R labels the unit cell of the WF, V is the volume of the unit cell, and U_mn^((k)) is an arbitrary unitary matrix. To construct maximally-localized WFs, U_mn^((k)) is chosen to minimize the following spread functional:
Ω=∑_n▒[〈r^2 〉_n-¯r_n^2 ] 										
where r ̅_n= ⟨0n│r│0n⟩ and 〈r^2 〉_n= ⟨0n│r^2│0n⟩. The minimization proceeds iteratively, based on an initial guess of localized orbitals.
For the case of interest in this work, where we wish to describe both the valence and conduction bands near the Fermi level, it is necessary to first select a set of bands to Wannierize, which includes separating the conduction bands from the free-electron-like bands that generally overlap with them in energy. The procedure to determine this localized subspace of Bloch wavefunctions proceeds similarly to minimization described above, where after an initial guess, the subspace is iteratively updated in order to minimize the spread function. After this initial disentanglement step, the Wannierization of the selected subspace proceeds as described above.
Due to the iterative non-linear minimization employed during both the disentanglement and Wannierization steps, the localization and utility of the final Wannier functions depend in practice on the initial choice of orbitals that are used to begin the disentanglement procedure, and which are then used as the initial guess for the Wannierization. Our initial guesses consist of a set of atomic orbitals we have chosen to describe all the chemically relevant orbitals for each element in typical elemental systems and compounds. We provide the list of the orbitals we select for each element in Table.
Table: Semi-core states used during Wannierization. Note that several lanthanides and actinides are avoided.
Element	POTCARs	N_electrons	Excluded_electrons	Projections	N_Wan_projections
Ag	Ag	11	0	s_d	6
Al	Al	3	0	s_p	4
Ar	Ar	8	0	s_p	4
As	As	5	0	s_p	4
Au	Au	11	0	s_d	6
B	B	3	0	s_p	4
Ba	Ba_sv	2	8	s_d	6
Be	Be_sv	2	2	s_p	4
Bi	Bi	5	0	s_p	4
Br	Br	7	0	s_p	4
C	C	4	0	s_p	4
Ca	Ca_sv	2	8	s_d	6
Cd	Cd	2	10	s_d	6
Ce	Ce	4	8	f_d_s	13
Cl	Cl	5	2	p	3
Co	Co	9	0	s_d	6
Cr	Cr_pv	6	6	s_d	6
Cs	Cs_sv	1	8	s_d	6
Cu	Cu_pv	11	6	s_d	6
Dy	Dy_3	12	3	s_f	8
Er	Er_3	14	5	f_s	8
Eu	Eu	9	8	f_s	8
F	F	5	2	p	3
Fe	Fe_pv	8	6	s_d	6
Ga	Ga_d	3	10	s_p	4
Gd	Gd	10	8	f_d_s	8
Ge	Ge_d	4	10	s_p	4
H	H	1	0	s	1
He	He	2	0	s	1
Hf	Hf_pv	4	6	s_d	6
Hg	Hg	12	0	s_p_d	9
I	I	7	0	s_p	4
In	In_d	3	10	s_p	4
Ir	Ir	9	0	s_d	6
K	K_sv	1	8	s_d	6
Kr	Kr	8	0	s_p	4
La	La	3	8	s_d_f	13
Li	Li_sv	1	2	s	1
Lu	Lu_3	17	8	f_d_s	13
Mg	Mg_pv	2	6	s_p	4
Mn	Mn_pv	7	6	s_d	6
Mo	Mo_pv	6	6	s_d	6
N	N	3	2	p	3
Na	Na_pv	1	6	s_p	4
Nb	Nb_pv	5	6	s_d	6
Nd	Nd_3	6	5	f_s	8
Ne	Ne	8	0	s_p	4
Ni	Ni_pv	10	6	s_d	6
O	O	4	2	p	3
Os	Os_pv	8	6	s_d	6
P	P	5	0	s_p	4
Pb	Pb_d	4	10	s_p	4
Pd	Pd	10	0	s_d	6
Pt	Pt	10	0	s_d	6
Rb	Rb_sv	1	8	s_d	6
Re	Re_pv	7	6	s_d	6
Rh	Rh_pv	9	6	s_d	6
Ru	Ru_pv	8	6	s_d	6
S	S	4	2	p	3
Sb	Sb	5	0	s_p	4
Sc	Sc_sv	3	8	s_d	6
Se	Se	6	0	s_p	4
Si	Si	4	0	s_p	4
Sm	Sm_3	8	3	f_s	8
Sn	Sn_d	4	10	s_p	4
Sr	Sr_sv	2	8	s_d	6
Ta	Ta_pv	5	6	s_d	6
Tb	Tb_3	9	0	f_s	8
Tc	Tc_pv	7	6	s_d	6
Te	Te	6	0	s_p	4
Th	Th	4	8	d_s	6
Ti	Ti_pv	4	6	s_d	6
Tl	Tl_d	3	10	s_p	4
U	U	9	5	f_s	8
V	V_pv	5	6	s_d	6
W	W_pv	6	0	s_d	6
Xe	Xe	8	0	s_p	4
Y	Y_sv	3	8	s_d	6
Zn	Zn	12	0	s_p_d	9
Zr	Zr_sv	4	8	s_d	6

 For many specific materials, it may be possible to select a smaller set of orbitals while still maintaining high-quality WFs that describe the bands of interest; however, our fairly inclusive set of orbitals is able Wannierize nearly all compounds in a high-throughput manner without human intervention. Because most applications of WFs are computationally inexpensive compared to the DFT calculations used to construct the WFs, in practice, our larger Wannier basis has only minimal computational cost. However, it is necessary to have enough empty bands in the underlying DFT calculation such that any empty orbitals chosen are included in the Bloch basis. We do not include any semicore orbitals in our Wannier basis, as they are generally well-separated in energy from the valence orbitals and are not necessary to describe bands near the Fermi level, and we exclude semicore bands from the disentanglement.
 During the disentanglement step, it is possible to choose an energy range that is included exactly (“the frozen window”), with the result that the Wannier band structure will exactly match the DFT band structure in this energy range and at the grid of k-points used in the Wannierization. We use a frozen window of ± 2 eV around the Fermi-energy. This window ensures that bands near the Fermi level are well described by the WTBH. For cases where the original WFs were unsatisfactory (see below), we found that lowering the lower bound of this window to include all of the valence bands often improves that WTBH, which we use as a second Wannierization setting. In order to validate our WTBH, we calculate the maximum absolute difference (μ) between the Wannier and DFT eigenvalues within an energy range of ± 2eV around the Fermi level:		
μ=(_nk^max) (|E_nk^DFT-E_nk^WTB |) 									
As discussed above, at the grid of k-points used in the construction of the WFs and within the frozen window, the eigenvalues should match exactly by construction. Therefore, to test the WTBH, we evaluate Eq. on the dense lines of k-points used to generate our band structure plots, which tests the WFs out of sample. A weakness of this evaluation method is that highly dispersive energy bands (high (dE_nk)/dk) can result in high μ values even if the WTBH is of good quality because any slight shift in the k-direction of a dispersive band will result in a large energy error. We consider that systems with μ less than 0.1 eV to useful for most applications, and we provide data for the user to evaluate individual WTBH for their own applications.
