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



###############
## Libraries ##
###############

from .__release import version
__version__ = version

## Abstract class
from .__DynamicalSystem import DynamicalSystem
from .__DiffDynSyst     import DiffDynSyst
from .__DiscDynSyst     import DiscDynSyst


## Continuous dynamical system
from .__Lorenz63 import Lorenz63
from .__Lorenz84 import Lorenz84
from .__Rossler  import Rossler


## Discrete dynamical system
from .__Henon import Henon
from .__Ikeda import Ikeda
from .__Mira  import Mira


## Functions
from .__dynamical_local_indexes import dynamical_local_indexes

