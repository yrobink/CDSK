
################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## This software is a computer program that is part of the CDSK (Chaotic      ##
## Dynamical System Kit) library. This library makes it possible              ##
## to generate classic (continuous and discrete) attractors, generate the     ##
## Mandelbrot and Julia set, and fit the local dimension.                     ##
##                                                                            ##
## This software is governed by the CeCILL-C license under French law and     ##
## abiding by the rules of distribution of free software.  You can  use,      ##
## modify and/ or redistribute the software under the terms of the CeCILL-C   ##
## license as circulated by CEA, CNRS and INRIA at the following URL          ##
## "http://www.cecill.info".                                                  ##
##                                                                            ##
## As a counterpart to the access to the source code and  rights to copy,     ##
## modify and redistribute granted by the license, users are provided only    ##
## with a limited warranty  and the software's author,  the holder of the     ##
## economic rights,  and the successive licensors  have only  limited         ##
## liability.                                                                 ##
##                                                                            ##
## In this respect, the user's attention is drawn to the risks associated     ##
## with loading,  using,  modifying and/or developing or reproducing the      ##
## software by the user in light of its specific status of free software,     ##
## that may mean  that it is complicated to manipulate,  and  that  also      ##
## therefore means  that it is reserved for developers  and  experienced      ##
## professionals having in-depth computer knowledge. Users are therefore      ##
## encouraged to load and test the software's suitability as regards their    ##
## requirements in conditions enabling the security of their systems and/or   ##
## data to be ensured and,  more generally, to use and operate it in the      ##
## same conditions as regards security.                                       ##
##                                                                            ##
## The fact that you are presently reading this means that you have had       ##
## knowledge of the CeCILL-C license and that you accept its terms.           ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
################################################################################
##                                                                            ##
## Copyright Yoann Robin, 2019                                                ##
##                                                                            ##
## yoann.robin.k@gmail.com                                                    ##
##                                                                            ##
## Ce logiciel est un programme informatique faisant partie de la librairie   ##
## CDSK (Chaotic Dynamical System Kit). Cette librairie permet de générer les ##
## attracteurs classiques (discret comme continue), de générer l'ensemble de  ##
## Julia et de Mandelbrot et d'estimer les dimensions locales.                ##
##                                                                            ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et  ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez    ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions     ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA   ##
## sur le site "http://www.cecill.info".                                      ##
##                                                                            ##
## En contrepartie de l'accessibilité au code source et des droits de copie,  ##
## de modification et de redistribution accordés par cette licence, il n'est  ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le     ##
## titulaire des droits patrimoniaux et les concédants successifs.            ##
##                                                                            ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques      ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au      ##
## développement et à la reproduction du logiciel par l'utilisateur étant     ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à      ##
## manipuler et qui le réserve donc à des développeurs et des professionnels  ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les    ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du     ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la       ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,   ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         ##
##                                                                            ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez     ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les  ##
## termes.                                                                    ##
##                                                                            ##
################################################################################
################################################################################

###############
## Libraries ##
###############


###############
## Functions ##
###############

## Ikeda {{{

#' Ikeda
#'
#' @description
#' Ikeda dynamical system
#'
#' @details
#' Discrete dynamical system.
#'
#' @examples
#' ike = CDSK::Ikeda$new( size = 200 )
#' X = ike$orbit(1000) 
#' ## X is an array with shape = (1000,200,2)
#' ## Each X[,i,] is an orbit
#' ## Each X[i,,] is a snapshot
#'
#' @export
Ikeda = R6::R6Class( "Ikeda" , 
	
	inherit = CDSK::DiscDynSystem,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	#' @field R  [float] Default = 1
	R = 1.,
	#' @field C1 [float] Default = 0.4
	C1 = 0.4,
	#' @field C2 [float] Default = 0.9
	C2 = 0.9,
	#' @field C3 [float] Default = 6.
	C3 = 6.,
	
	#################
	## Constructor ##
	#################
	
	#' @description
    #' Create a new Ikeda object.
	#' @param R  [float] Default = 1
	#' @param C1 [float] Default = 0.4
	#' @param C2 [float] Default = 0.9
	#' @param C3 [float] Default = 6.
	#' @param size [integer] Number of initial condition simultaneously solved
	#' @return A new `Ikeda` object.
	initialize = function( R = 1. , C1 = 0.4 , C2 = 0.9 , C3 = 6. , size = 1 )
	{
		super$initialize( 2 , size , base::matrix( base::c( 0 , -1 , 1 , 0 ) , nrow = 2 , ncol = 2 ) )
		self$R  = R
		self$C1 = C1
		self$C2 = C2
		self$C3 = C3
	}
	
	
	#############
	## Methods ##
	#############
	
	),
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	#############
	## Methods ##
	#############
	
	equation = function( t , X , par = NULL )
	{
		Xnext = numeric(length(X))
		tau = self$C1 - self$C3 / ( 1. + X[private$i[[1]]]^2 + X[private$i[[2]]]^2 )
		Xnext[private$i[[1]]] = self$R + self$C2 * ( X[private$i[[1]]] * base::cos(tau) - X[private$i[[2]]] * base::sin(tau) )
		Xnext[private$i[[2]]] = self$C2 * ( X[private$i[[2]]] * base::cos(tau) + X[private$i[[1]]] * base::sin(tau) )
		return(Xnext)
	}
	
	)
)
##}}}

