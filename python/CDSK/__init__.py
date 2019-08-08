# -*- coding: utf-8 -*-

###############
## Libraries ##
###############

__version__ = "0.1.0a0"

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
from .__LocalDimension import localDimension


