
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

library(devtools)
library(roxygen2)


## Read command line arguments
##============================

args = commandArgs( trailingOnly = TRUE )
verbose = FALSE
install = FALSE
build   = TRUE
check   = FALSE
if( length(args) > 0 )
{
	for( a in args )
	{
		if( a == "-v" || a == "--verbose" )
		{
			verbose = TRUE
		}
		if( a == "-i" || a == "--install" )
		{
			install = TRUE
		}
		if( a == "-nb" || a == "--not-build" )
		{
			build = FALSE
		}
		if( a == "-c" || a == "--check" )
		{
			check = TRUE
		}
	}
}


## Building
##=========
cdsk = ""
if( build )
{
	if( verbose ) cat( "Generation of Rd files with roxygen" )
	roxygen2::roxygenize("CDSK")
	if( verbose ) cat( "Load of CDSK to generate Rd files with Rcpp" )
	devtools::load_all("CDSK")
	if( verbose ) cat( "Generation of Rd files for cpp with roxygen" )
	roxygen2::roxygenize("CDSK")
	if( verbose ) cat( "Final build" )
	cdsk = devtools::build("CDSK")
}


## Check
##======
if( check )
{
	if( verbose ) cat( "Check CDSK" )
	devtools::check( "CDSK" )
}


## Installation
##=============

if( install )
{
	if( verbose ) cat( "Installation" )
	if( cdsk == "" )
	{
		files = base::list.files()
		cdsk = ""
		
		for( f in files )
		{
			f_split = base::unlist(base::strsplit( f , "[.]" ))
			if( length(f_split) > 2 )
			{
				f_split = base::rev(f_split)
				if( f_split[1] == "gz" && f_split[2] == "tar" )
				{
					cdsk = f
				}
			}
		}
	}
	if( cdsk == "" )
	{
		cat( "Error, CDSK not built, so can not be installed" )
	}
	else
	{
		install.packages(cdsk)
	}
}

