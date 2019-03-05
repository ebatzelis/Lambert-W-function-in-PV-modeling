# Lambert-W-function-in-PV-modeling
This repository contains the MATLAB code and C code for 6 implementation approaches of the Lambert W function in photovoltaic (PV) models, discussed in the paper:
Efstratios Batzelis, Georgios Anagnostou, Chandan Chakraborty and Bikash Pal, “Computation of the Lambert W function in photovoltaic modelling,” accepted for presentation in ELECTRIMACS 2019.
The MATLAB code includes .m files for: the Asymptotic (7 and 4terms), the Hybrid, the Simple and the Analytical formulae. The C code codes contains the aforementioned functions plus the iterative Halley’s method adopted in MATLAB lambertw function.
For applications to microcontrollers, please use the C code, as it is appropriately optimized (minimum number of divisions, logarithmic/exponential evaluations, factoring terms etc.).
Please share and cite the paper above.
For any queries: e.batzelis@imperial.ac.uk

