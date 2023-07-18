Introduction
---
The code is a parallel Fortran program targeted at the 1D forward modeling specifically designed for overhead transmission lines.
It can accurately calculate electromagnetic fields excited by overhead transmission lines under various factors, including the design parameters, sag, and orientation of transmission lines.

Prerequisites
---
Before using this code, you need to install **gfortran compiler** and **OpenMP**.

Files
---
The program consists of five **f90** files, which are as follows:
- ***main.f90***: The main program, which is used for data preparation, module invocation, and result output.
- ***module_coor_trans_intrans.f90***: This file is used for coordinate translation and inverse coordinate translation.
- ***module_freq_forward.f90***: This file culculates forward responses.
- ***module_GS_HK.f90***: This file saves Hankel sampling points and weights, and calculates Gauss points and Gauss weights.
- ***module_parameter_set.f90***: This file sets parameters used during forward modeling.

Usage
---
To run the program, the user needs to provide three files: the model file, the data template file, and the configuration file. Then use the flowing command to run emfem:

Author
---
- Name: Siwei Zhu
- Affiliation: China University of Geosciences, Wuhan, China
- Email: <selwey1996@gmail.com>

License
---
This code is distributed under GPL License. Please see the file [LICENSE](./LICENSE) for more details.
