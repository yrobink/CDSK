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

## Read release informations
from CDSK.__release import name
from CDSK.__release import version
from CDSK.__release import description
from CDSK.__release import long_description
from CDSK.__release import authors
from CDSK.__release import authors_email
from CDSK.__release import src_url
from CDSK.__release import license


## Required elements
author           = ", ".join(authors)
author_email     = ", ".join(authors_email)
packages         = ["CDSK","CDSK.fractal"]
package_dir      = { "CDSK" : "CDSK" }
requires         = [ "numpy",
					 "scipy",
					 "matplotlib",
					 "sklearn",
					 "SDFC (>=0.5.0)"]
keywords         = ["Dynamical-system","chaos","attractor","local-dimension"]
platforms        = ["linux","macosx"]

## Now the setup
from distutils.core import setup

setup(  name             = name,
		version          = version,
		description      = description,
		long_description = long_description,
		author           = author,
		author_email     = author_email,
		url              = src_url,
		packages         = packages,
		package_dir      = package_dir,
		requires         = requires,
		license          = license,
		keywords         = keywords,
		platforms        = platforms
     )
