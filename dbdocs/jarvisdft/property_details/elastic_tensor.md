We use finite difference method on conventional cells of systems for ET calculations. For bulk material, the compliance tensor can be obtained by:
〖s_ij=C〗_ij^(-1)											     
Now, several other elastic properties calculated from Cij and sij. Some of the important properties are given below:
KV = ((C11+C22+C33) + 2(C12+C23+C31))/9							      
GV = ((C11+C22+C33) − (C12+C23 + C31) +3 (C44+C55+C66))/15				       
KR = ((s11+s22+s33) + 2(s12+s23+s31))-1								      
GR = 15(4(s11+s22+s33) - 4(s12+s23+s31) + 3(s44+s55+s66))-1					      
KVRH =(KV+KR)/2										     
GVRH =(GV+GR)/2										     
ν = (3KVRH − 2GVRH)/(6KVRH+2GVRH))						                 
Here KV and GV  are Voigt bulk and shear modulus, and KR and GR Reuss-bulk and shear modulus respectively. The homogenous Poisson ratio is calculated as ν. The EC data can be also used to predict the ductile and brittle nature of materials with Pugh (Gv/Kv) and Pettifor criteria (C12-C44) . Materials with Pugh ratio value >0.571 and Pettifor criteria <0 should be brittle, while materials with Pugh ratio value <0.571 and Pettifor criteria >0 should be ductile.
For monolayer material calculations, the elastic tensor obtained from DFT code such as VASP, assumes periodic-boundary-condition (PBC). Therefore, cell vectors are used to calculate the area which again is used in computing stress. When dealing with the monolayer, an arbitrary vacuum padding is added in one of the direction (say z-direction). When computing EC we need to correct the output by eliminating the arbitrariness of the vacuum padding. We do that as a post-processing step by multiplying the Cij components  (i,j≠3) by the length of the vacuum padding. Therefore, the units of EC turn into Nm-1 from Nm-2. For example, in order to calculate C11 (stress computed in x direction), the area is computed using normal of y and z-vectors. Obviously, the z-vector is arbitrary, so if we multiply the output by z-vector magnitude we get rid of the arbitrariness of z and also get C11 in Nm-1. As shown in Fig. 1, the z-vector magnitude is the z-simulation. The above discussion can also be expressed as the following:
σ_VASP=F/A=├ F/|z||l| ┤|_(l∈(x,y) )									     
σ_mono=〖|z|× σ〗_VASP	
