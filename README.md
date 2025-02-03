# UltimateKalman: Flexible Kalman Filtering and Smoothing Using Orthogonal Transformations

This repository contains the source code of several Kalman filters and smoothers (including the repository's namesake, UltimateKalman) in three different programming languages: MATLAB, C and Java.

- Release v1.2.0, contains the sequential UltimateKalman algorithm, which is documented carefully in an article in the [ACM Transactions on Mathematical Software](https://doi.org/10.1145/3699958). The algorithm is implemented monolithically in all 3 languages. Please cite this article when citing the sequential UltimateKalman algorithm or its implementation. This release is also available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14543666.svg)](https://doi.org/10.5281/zenodo.14543666).  

- Release v2.0.0 contains two parallel Kalman smoothers, as well as a conventional sequential Kalman filter and smoother. These implementations are described in an article by Shahaf Gargir and Sivan Toledo accepted to IPDPS 2025. The actual parallel implementations are available only in C/C++, but there are also sequential Matlab implementations of all the algorithms, for clarity. The Matlab implementation is now modular and uses a class hierarchy, to eliminate code duplication (so it is no longer monolithic).

## User Guide

A detailed user guide (for version 1.2.0!) is available in [userManual/userManual.pdf](userManual/userManual.pdf).

## License

Copyright 2020-2024 Sivan Toledo.
 
 UltimateKalman is free software; you can redistribute it and/or modify
    it under the terms of either:

 the GNU Lesser General Public License as published by the Free
        Software Foundation; either version 3 of the License, or (at your
        option) any later version.

or

the GNU General Public License as published by the Free Software
        Foundation; either version 2 of the License, or (at your option) any
        later version.

or both in parallel, as here, or 

the Apache License as published by the Apache Software
        Foundation; either version 2.0 of the License, or (at your option) any
        later version.
        
WITH THE ADDITIONAL REQUIREMENT 
    that if you use this software or derivatives of it, directly or indirectly, to produce
    research that is described in a research paper, you need to cite the most
    up-to-date version of the article that describes UltimateKalman in your paper.
    
Currently, the version to cite is [the version on arXiv](https://arxiv.org/abs/2207.13526).

UltimateKalman is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

You should have received copies of the GNU General Public License and the
    GNU Lesser General Public License along with this software.  If not,
    see https://www.gnu.org/licenses/.
    
