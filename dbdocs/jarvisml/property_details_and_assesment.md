Currently, there are two types of data that are machine-learned in JARVIS-ML: discrete and image-based. The discrete target is obtained from the JARVIS-DFT database for 3D and 2D materials. There have been several descriptor developments as attempts to capture the complex chemical-structural information of a material. We compute CFID descriptors for most crystal structures in various databases. Many of these structures are non-unique but can still be used for pre-screening applications. The CFID can also be applied to other materials classes such as molecules, proteins, point defects, free surfaces, and heterostructures, which are currently ongoing projects. These descriptor datasets, along with JARVIS-DFT and other databases, act as input and outputs for machine learning algorithms. The CFID consists of 1557 descriptors for each material: 438 average chemical, 4 simulation-box-size, 378 radial charge-distribution, 100 radial distribution, 179 angle-distribution up to first neighbor, and another 179 for the second neighbor, 179 dihedral angle up to fist neighbor and 100 nearest neighbor descriptors. We have converted at least 1 million atomic structures to CFID descriptors and that dataset will be made available soon. After generating the features, we trained several property models. 
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
