# CDSK (Chaotic Dynamical System Kit)


## Features
- Generation of classic chaotic attractors
- Generation of Mandelbrot and Julia set
- Estimation of Local Dimension and persistence

## Installation

Requires the [pmetric](https://github.com/yrobink/pmetric) packages.

The file *build.R* can be used to compile and install CDSK. For almost all users, it is sufficient to launch:

~~~~
Rscript build.R -v -c -i
~~~~

Available options are:
- *-v* or *--verbose*
- *-c* or *--check* to launch a check of CDSK before installation
- *-i* or *--install* to install
- *-nb* or *--not-build*, for developpers


## License

Copyright Yoann Robin, 2019, [CeCILL-C](https://cecill.info/licences.en.html) license (open source)

This software is a computer program that is part of the CDSK (Chaotic
Dynamical System Kit) library. This library makes it possible
to generate classic (continuous and discrete) attractors, generate the
Mandelbrot and Julia set, and fit the local dimension.

This software is governed by the CeCILL-C license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-C 
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-C license and that you accept its terms.

