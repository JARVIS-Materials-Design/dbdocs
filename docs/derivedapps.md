# JARVIS-OtherApps

## Overview

Based on the main JARVIS-apps such as JARVIS-DFT, we have developed several derived apps such as JARVIS-WannierTB, and JARVIS-Heterostructure etc.
These apps allow to go beyond the existing data in JARVIS and let users obtain predictions on their custom inputs. These apps are based Django
and the foundation codes to build these apps are provided at JARVIS-Tools GitHub page.

## JARVIS-WannierTB

 - JARVIS-WannierTB provides Wannier Tight binding Hamiltonians and an interface to solve these Hamiltonians for arbitrary k-points on-the-fly. 
 - The JARVIS-Wannier tight-binding Hamiltonians were generated for both 3D and 2D materials.
 - We evaluate the accuracy of the WTBHs by comparing the Wannier band structures to directly 
calculated DFT band structures on both the set of k-points used in 
the Wannierization as well as independent k-points from high symmetry lines.


## JARVIS-Heterostructure

 - JARVIS-Heterostructure provide interface-design app and tools for 2D materials in the JARVIS-DFT database using Zur algorithm and Anderson rule.
 - Some of the properties available are: work function, band-alignment, and heterostructure classification. 
 - The band alignment data can be also used for identifying photocatalysts and high-work function 2D-metals for contacts. 
 - We validate our results by comparing them to experimental data as well as hybrid-functional predictions.


## References
1. 	[Database of Wannier Tight-binding Hamiltonians using High-throughput Density Functional Theory, arXiv:2007.01205](https://arxiv.org/abs/2007.01205)
2. 	[Efficient Computational Design of 2D van der Waals Heterostructures: Band-Alignment, Lattice-Mismatch, Web-app Generation and Machine-learning, arXiv:2004.03025](https://arxiv.org/abs/2004.03025).

