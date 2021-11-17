# CDSK (Chaotic Dynamical System Kit)


## Features
- Generation of classic chaotic attractors
- Generation of Mandelbrot and Julia set
- Estimation of Local Dimension and persistence

## Installation

Require the [pmetric](https://github.com/yrobink/pmetric) package.

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

Copyright(c) 2021 Yoann Robin

CDSK is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CDSK is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with CDSK.  If not, see <https://www.gnu.org/licenses/>.

