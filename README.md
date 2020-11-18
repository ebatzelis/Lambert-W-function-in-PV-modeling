# Lambert-W-function-in-PV-modeling
This repository contains the MATLAB code and C code for 6 implementation approaches of the Lambert W function in photovoltaic (PV) models, discussed in the paper:

Efstratios Batzelis, Georgios Anagnostou, Chandan Chakraborty and Bikash Pal, “Computation of the Lambert W function in photovoltaic modelling,” in ELECTRIMACS 2019, Springer, 2020, pp. 583–595.
https://link.springer.com/chapter/10.1007%2F978-3-030-37161-6_44

The MATLAB code includes .m files for: the Asymptotic (7 and 4terms), the Hybrid, the Simple and the Analytical formulae. The C code contains the aforementioned functions plus the iterative Halley’s method adopted in the MATLAB lambertw function.

For implementation to a microcontroller, please use the C code as it is appropriately optimized (minimum number of divisions, logarithmic/exponential evaluations, factoring terms etc.).

Please share and cite the paper above. For any queries: stratis.batzelis@imperial.ac.uk

