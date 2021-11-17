# -*- coding: utf-8 -*-

## Copyright(c) 2021 Yoann Robin
## 
## This file is part of CDSK.
## 
## CDSK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## CDSK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with CDSK.  If not, see <https://www.gnu.org/licenses/>.

version_major = 1
version_minor = 0
version_patch = 0
version_extra = "a2"
version       = "{}.{}.{}{}".format(version_major,version_minor,version_patch,version_extra)


##

name = "CDSK"


##

description      = "Chaotic Dynamical System Kit"
long_description = \
"""
Chaotic Dynamical System Kit (CDSK)

- Generation of classic chaotic attractors
- Generation of Mandelbrot and Julia set
- Estimation of Local Dimension and persistence
"""


##

license = "GNU-GPL3"


##

authors       = ["Yoann Robin"]
authors_email = ["yoann.robin.k@gmail.com"]
author_doc    = ", ".join( [ ath + " ({})".format(athm) for ath,athm in zip(authors,authors_email) ] )


##

src_url = "https://github.com/yrobink/CDSK"



