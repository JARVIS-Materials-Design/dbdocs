The transport properties were calculated using the Boltzmann transport equation (BTE) implemented in the BoltzTrap code. The BTE is used to investigate the non-equilibrium behavior of electrons and holes by statistically averaging all possible quantum states given by the equation below:
df(k,T,t)/dt=(∂f(k,T,t)/∂t)_s-dk/dt ∇_k f(k,T,t)-v(k) ∇_r f(k,T,t).					
f is the electron distribution, which is a function of state k, temperature T and time t, and v(k) are the electron group velocities. The three terms on the right-hand side of Eq. refer, respectively, to the temporal rate of change of f due to all scattering processes, the rate of change of f due to external forces, and the diffusion from the carrier density gradient. If the external forces consist only of a low electric field, E, and no magnetic field, B, such that dk/dt=eE/h then the low-filed BTE is given by:
df(k,T,t)/dt+v(k) ∇_r f(k,T)+eE/h ∇_k f(k,T)=(∂f(k,T,t)/∂t)_s.				
Now, f can be described as a first-order (linear) perturbation from the (equilibrium) Fermi-Dirac distribution, f0, due to scattering
(∂f(k,T,t)/∂t)_s=-(f(k)-f_0 (k))/τ								
where 
f_0 [E(k)]=1/(e^([E(k)-E_F ]/k_B T)+1)							
where the dependence of E on k is given by the electronic band structure, and the various
scattering terms and time dependence are lumped into the electronic relaxation time, τ. The computation of the relaxation time is very computationally expensive, especially in a high-throughput context.
