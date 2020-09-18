 We carry out the Density functional perturbation theory (DFPT) (as implemented in the VASP code) calculation on the standard conventional cell for each material. We determine the Born-effective charge, piezoelectric and dielectric (ionic+electronic parts) tensors and the phonon eigenvectors. DFPT is also used to calculate Infrared and Raman intensities for selected materials. As mentioned earlier, it is important to converge K-points, cut-off and parameters to get reliable results which is taken into account in the JARVIS-DFT. Given an insulating system with N atoms per cell, with cell volume Ω0, atomic displacements um (m=(1…3N)), homogenous strain ηj (j=(1…6)), homogenous electric fields  Eα (α=(x,y,z)), energy E, the force-constant matrix (K_mn), internal strain tensor (Λ_mj),  the  dielectric susceptibility (χ_αβ), Born dynamical effective charge (Z_mα), piezoelectric stress tensor (e_αj) and piezoelectric strain tensor (d_αj) are calculated as follows (SI units are used throughout)65-69: 
K_mn=Ω_0  (∂^2 E)/(∂u_m ∂u_n ) |_(E,η)									
Z_mα=-Ω_0  (∂^2 E)/(∂u_m ∂E_α ) |_η									
Λ_mj=-Ω_0  (∂^2 E)/(∂u_m ∂η_j ) |_E								
χ_αβ==-(∂^2 E)/(∂E_α ∂E_β ) |_(u,η)+Ω_0^(-1) Z_mα (K^(-1) )_mn Z_nβ					
e_αj=-(∂^2 E)/(∂E_α ∂η_j ) |_u+Ω_0^(-1) Z_mα (K^(-1) )_mn Λ_nj						
d_αj=S_jk^((E) ) e_αk									
The dielectric constant can be derived from the dielectric susceptibility using:
〖ε_αβ=ε_0 (δ_αβ+χ〗_αβ)									
The first term represents the electronic contribution and the second term the ionic contribution for DL and PZ constants respectively.
The PZ is a 3x6 tensor, the DL 3x3 and the BEC Nx3x3 tensor. The IR intensity of phonon modes is calculated using:
f(n)=∑_α▒|∑_sβ▒Z_αβ  (s) e_β (s,n)|^2 							
where e_β (s,n) is the normalized vibrational eigenvector of the nth phonon mode of the sth atom in the unit cell, and α, β are the cartesian coordinates. Z_αβ (s) is the Born effective charge tensor of sth atom (here we explicitly write both the cartesian indices of Z). These approaches are universal and have been already applied to various material classes. More details about the DFPT formalism can be found in elsewhere. 
