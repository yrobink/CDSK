
###############
## Libraries ##
###############


###############
## Functions ##
###############


## Julia

#' Julia
#'
#' Julia set
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param c [complex number]
#'        Constant in the iteration of Julia set
#' @param x [vector]
#'        x-coordinates of Julia set
#' @param y [vector]
#'        y-coordinates of Julia set
#' @param classic_set [string]
#'        A pre-defined set, values are "set0", "set1" , "set2" , "set3" or "set4".
#' @param nx [integer]
#'        Numbers of step in x-axis
#' @param ny [integer]
#'        Numbers of step in y-axis
#' @param maxit [integer]
#'        Number of iterations to estimate the Julia set
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(c,x,y,classic_set,nx,ny,maxit)}}{This method is used to create object of this class with \code{Julia}}
#'   \item{\code{run()}}{Estimate the ratio of the Julia set}
#'   \item{\code{plot()}}{Plot}
#' }
#' @examples
#' m = Julia$new( classic_set = "set0" , nx = 200 , ny = 200 )
#' m$run()
#'
#' @export
Julia = R6::R6Class( "Julia" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	ratio = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( c = NULL , x = NULL , y = NULL , classic_set = NULL , nx = NULL , ny = NULL , maxit = 200 )
	{
		private$c = c
		private$x = x
		private$y = y
		private$nx = nx
		private$ny = ny
		private$maxit = maxit
		
		if( !is.null(classic_set) && ( !is.null(nx) || !is.null(ny) ) )
		{
			private$nx = if( is.null(nx) ) ny else nx
			private$ny = if( is.null(ny) ) nx else ny
			
			if(      classic_set == "set0" ) private$set_set0()
			else if( classic_set == "set1" ) private$set_set1()
			else if( classic_set == "set2" ) private$set_set2()
			else if( classic_set == "set3" ) private$set_set3()
			else if( classic_set == "set4" ) private$set_set4()
		}
		else if( !is.null(x) && !is.null(y) )
		{
			private$nx = length(x)
			private$ny = length(y)
		}
		else
		{
			print("Error, Need (x,y) set or (classic_set,nx [,ny]) set!" )
		}
	},
	
	
	#############
	## Methods ##
	#############
	
	run = function()
	{
		grid = as.matrix( base::expand.grid( private$y , private$x ) )
		z = complex( 1 , grid[,2] , grid[,1] )
		ratio = numeric( private$nx * private$ny )
		
		for( i in 1:private$maxit )
		{
			z = z^2 + private$c
			idx = base::abs(z) < 2
			ratio[idx] = ratio[idx] + 1
		}
		
		self$ratio = matrix( ratio , nrow = private$nx , ncol = private$ny ) / private$maxit
	},
	
	plot = function()
	{
		grDevices::dev.new()
		graphics::image( self$ratio , col = base::rev(heat.colors(12)) )
	}
	
	),
	
	
	private = list(
	
	###############
	## arguments ##
	###############
	
	c = NULL,
	x = NULL,
	y = NULL,
	nx    = NULL,
	ny    = NULL,
	maxit = NULL,
	
	
	#############
	## Methods ##
	#############
	
	set_set0 = function()
	{
		private$c = complex(1,0.3,0.5)
		private$x = base::seq( 0.25 , 0.53 , length = private$nx )
		private$y = base::seq( 0.4  , 0.72  , length = private$ny )
	},
	
	set_set1 = function()
	{
		private$c = complex(1,0.285,0.01)
		private$x = base::seq( -1 , 1 , length = private$nx )
		private$y = base::seq( -1.1  , 1.1  , length = private$ny )
	},
	
	set_set2 = function()
	{
		private$c = complex(1,0.285,0.013)
		private$x = base::seq( -1 , 1 , length = private$nx )
		private$y = base::seq( -1.1  , 1.1  , length = private$ny )
	},
	
	set_set3 = function()
	{
		private$c = complex( 1 , -0.038088 , 0.9754633 )
		private$x = base::seq( -2e-3 , 2e-3 , length = private$nx )
		private$y = base::seq( -2e-3 , 2e-3 , length = private$ny )
	},
	
	set_set4 = function()
	{
		private$c = complex( 1 , -0.8 , 0.156 )
		private$x = base::seq( -1.55 , 1.55 , length = private$nx )
		private$y = base::seq( -1    , 1    , length = private$ny )
	}
	
	)
)

