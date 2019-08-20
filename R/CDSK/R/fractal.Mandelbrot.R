
###############
## Libraries ##
###############


###############
## Functions ##
###############


## Mandelbrot

#' Mandelbrot
#'
#' Mandelbrot set
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'        x-coordinates of Mandelbrot set
#' @param y [vector]
#'        y-coordinates of Mandelbrot set
#' @param classic_set [string]
#'        A pre-defined set, values are "set0", "set1" , "set2" , "set3" or "set4".
#' @param nx [integer]
#'        Numbers of step in x-axis
#' @param ny [integer]
#'        Numbers of step in y-axis
#' @param maxit [integer]
#'        Number of iterations to estimate the Mandelbrot set
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(x,y,classic_set,nx,ny,maxit)}}{This method is used to create object of this class with \code{Mandelbrot}}
#'   \item{\code{run()}}{Estimate the ratio of the Mandelbrot set}
#'   \item{\code{plot()}}{Plot}
#' }
#' @examples
#' m = Mandelbrot$new( classic_set = "set0" , nx = 200 , ny = 200 )
#' m$run()
#'
#' @export
Mandelbrot = R6::R6Class( "Mandelbrot" ,
	
	public = list(
	
	###############
	## arguments ##
	###############
	
	ratio = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( x = NULL , y = NULL , classic_set = NULL , nx = NULL , ny = NULL , maxit = 200 )
	{
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
		c = complex( 1 , grid[,2] , grid[,1] )
		z = numeric( private$nx * private$ny )
		ratio = numeric( private$nx * private$ny )
		
		for( i in 1:private$maxit )
		{
			z = z^2 + c
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
		private$x = base::seq( -0.62 , -0.42 , length = private$nx )
		private$y = base::seq( -0.7  , -0.5  , length = private$ny )
	},
	
	set_set1 = function()
	{
		private$x = base::seq( -0.6 , -0.595 , length = private$nx )
		private$y = base::seq( -0.665  , -0.66  , length = private$ny )
	},
	
	set_set2 = function()
	{
		private$x = base::seq( -0.57 , -0.56 , length = private$nx )
		private$y = base::seq( -0.646  , -0.636  , length = private$ny )
	},
	
	set_set3 = function()
	{
		private$x = base::seq( -0.5663 , -0.5656 , length = private$nx )
		private$y = base::seq( -0.6394  , -0.6387  , length = private$ny )
	},
	
	set_set4 = function()
	{
		private$x = base::seq( -1.5 , -0.5 , length = private$nx )
		private$y = base::seq( -0.5  , 0.5  , length = private$ny )
	}
	
	)
)

