
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include "Pairwise_distances.hpp"


RCPP_MODULE(yada){
	Rcpp::function( "cpp_pairwise_distances_XYstr"  , &cpp_pairwise_distances_XYstr , "" ) ;
	Rcpp::function( "cpp_pairwise_distances_Xstr"   , &cpp_pairwise_distances_Xstr , "" ) ;
	Rcpp::function( "cpp_pairwise_distances_XYCall" , &cpp_pairwise_distances_XYCall , "" ) ;
	Rcpp::function( "cpp_pairwise_distances_XCall"  , &cpp_pairwise_distances_XCall , "" ) ;
}

