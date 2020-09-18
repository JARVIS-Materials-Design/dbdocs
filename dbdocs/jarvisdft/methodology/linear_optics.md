To obtain the optical properties of the materials, we calculated the imaginary part of the dielectric function from the Bloch wavefunctions and eigenvalues (neglecting local field effects). We introduced three times as many empty conduction bands as valance bands. This treatment is necessary to facilitate proper electronic transitions. We choose 5000 energy grid points to have a sufficiently high resolution in dielectric function spectra. The imaginary part is calculated as:
ε_αβ^((2)) (E)=(4π^2 e^2)/Ω^2    lim┬(q→0)⁡〖1/q^2  ∑_(c,v,k ⃗)▒〖2w_k ⃗  δ(ξ_(ck ⃗ )-ξ_(vk ⃗ )-E)〈Ψ_(ck ⃗+e ⃗_α q) |Ψ_(vk ⃗ ) 〉 〈Ψ_(vk ⃗ ) |Ψ_(ck ⃗+e ⃗_β q) 〉^* 〗〗	
where e is electron charge, Ω is the cell volume, E the energy,   is the Fermi-weight of each k-point,  are unit vectors along the three Cartesian directions,   is the cell-periodic part of the pseudo-wavefunction for band n and k-point k, q stands for the Bloch vector of an incident wave, c and v stand for conduction and valence bands,   stands for eigenvalues of the corresponding bands respectively. The matrix elements on the right side of Eq. (1) capture the transitions allowed by symmetry and selection rules. The real part of the dielectric tensor ε_αβ^((1)) is obtained by the usual Kramers-Kronig transformation:
ε_αβ^((1)) (E)=1+2/π P∫_0^∞▒(ε_αβ^((2)) (E^' ) E^')/((E^' )^2-E^2+iη) dE^'			     				  
where P denotes the principle value, and η is the complex shift parameter taken as 0.1. Moreover, as the dielectric function is a tensorial quantity, we use the crystallographic average of the dielectric function (written as ε^((1))  and ε^((2))), obtained by diagonalizing the dielectric tensor for each energy and averaging the diagonal elements.
Using, ε^((1)) and ε^((2)) the absorption coefficient α(E) is defined as:
α(E)=2E/ℏc √((√((ε^((1)) (E))^2+(ε^((2)) (E))^2 )-(ε^((1)) (E)))/2)							   
where c is the speed of light.
Next, the SLME (ɳ) is defined as the ratio of the maximum output power density (Pmax) and the total incident solar energy density (Pin). Pmax is obtained by numerically maximizing the product of current density J and voltage V.
ɳ=P_max/P_in  										
Assuming the solar cell at temperature T behaves as an ideal diode and is illuminated under the photon flux Isun, J and V follow the following equation:
J=J_sc-J_0 (e^(eV/kT)-1)								
where e is the elementary charge, V the potential over the absorber layer and k is Boltzmann’s constant. The first term is the short-circuit current density Jsc given by:
J_sc=e∫_0^∞▒〖a(E) I_sun (E)dE〗									
where and a(E) is the photon absorptivity, I_sun  is the AM1.5G solar spectrum. The a(E) depends on the absorption coefficient (α) (Eq. (3)) and thickness (L) of the material.
a(E)=1-e^(-2α(E)L)										
The coefficient of the second term in Eq. (5) is the reverse saturation current (J0), which corresponds to the total (radiative and non-radiative) electron-hole recombination current at equilibrium in the dark:
J_0^ =J_0^r+J_0^nr=(J_0^r)/f_r 										
Here, fr is defined as the fraction of the radiative recombination current. For the SLME, fr is approximated using:
f_r=e^(((E_g-E_g^da)/kT) )											
Where E_g is the fundamental and E_g^da  is the direct allowed bandgap of a material.
Following the principle of detailed balance, the rates of emission and absorption through cell surfaces must be equal in equilibrium in the dark. Hence, J_0^r can be calculated from the rate at which black-body photons from the surrounding thermal bath are absorbed through the front surface, given by:
J_0^r=eπ∫_0^∞▒〖a(E) I_bb (E,T)dE〗									
where Ibb is the black-body spectrum at temperature T. Both the solar spectrum I_sun  and black-body spectrum I_bb are expressed in terms of the photon flux.
In order to maximize the power density, Eq.  can be re-written as:
ɳ=P_max/P_in =max⁡〖{(J_sc-J_0 (e^(eV/kT)-1))V}_V 〗/(∫_0^∞▒〖EI_sun (E)dE〗)								
Therefore, the material-property related inputs in calculating the SLME are α(E),f_r,L and T. In this work, we assume material thickness (L) as 500 nm and temperature (T) as 300 K. 
VASP uses a complex shift (CSHIFT) in the Kramers-Kronig relation to calculate the real part of the dielectric tensor, and also determines the corresponding imaginary part for consistency.
