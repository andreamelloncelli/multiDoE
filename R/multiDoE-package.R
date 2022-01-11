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
#' Although the implemented methods are designed to handle designs with
#' experimental factors in at least two different strata, their flexibility
#' allows their application to even the simplest cases of completely randomized
#' and randomized block designs. Whatever the case of interest, the designs
#' manageable by the package are balanced. The user can choose the structure of
#' the experimental design by defining:
#' \itemize{
#' \item the number of strata of the experiment;
#' \item the number of experimental factors in each stratum;
#' \item the number of experimental units in each stratum and thus the total number of runs;
#' \item the number of levels of each experimental factor;
#' \item the presence or not of blocking factors.
#' }
#' It is possible to choose the a priori model among: the model with main effects
#' only, the model with main and interaction effects and the full quadratic model.
#' Finally, estimation of the ratios of error variances in consecutive strata is
#' required. With the package \code{MultiDoE} it is possible to obtain experimental
#' designs that optimize any combination of the following criteria: "I", "Id",
#' "D", "Ds", "A" and "As". Depending on the function used, it is possible to obtain
#' either a single optimal solution for the optimization problem of interest or
#' the set of solutions belonging to the (approximate) Pareto front. Once the
#' Pareto front is available, the package offers: a method for the graphical
#' visualization of the Pareto front (up to a three-dimensional criteria space);
#' an objective method for selecting the best design with respect to the entire
#' set of criteria considered and an objective method for selecting the best design
#' with respect to to the criteria considered individually.
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
