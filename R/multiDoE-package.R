# Roxygen docs:
# https://roxygen2.r-lib.org/articles/rd.html#packages

#' @description The R package \code{MultiDoE} can be used to construct
#' multi-stratum experimental designs that optimize up to six statistical
#' criteria simultaneously. The algorithms implemented in the package to solve
#' such optimization problems are the \emph{MS-Opt} and \emph{MS-TPLS} algorithms
#' proposed in Sambo, Borrotti, Mylona, and Gilmour (2017). The former relies on
#' a search procedure to find locally optimal solutions; the latter, by embedding
#' MS-Opt in a TPLS framework, is able to generate a good Pareto front
#' approximation for the optimization problem under study.
#'
#'
#' @details
#' XXX TODO FIX ME XXXX The only function you're likely to need from roxygen2 is [roxygenize()].
#' Otherwise refer to the vignettes to see how to format the documentation.
#' This comes from the R/multiDOE-package.R
#'
#' @references
#' M. Borrotti and F. Sambo and K. Mylona and S. Gilmour. A multi-objective
#' coordinate-exchange two-phase local search algorithm for multi-stratum
#' experiments. Statistics & Computing, 2017.
#'
"_PACKAGE"
