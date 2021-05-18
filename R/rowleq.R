# check if row A is less than or equal to row B

# function res = rowleq( A, B )
# i = 1;
# while (A(i) == B(i)) && (i < length(A))
# i = i + 1;
# end
#
# res = (A(i) <= B(i));


rowleq <- function(A, B) {
  argmin = which.min(A==B)   # restituisce il primo FALSE (0)
  return(A[argmin] <= B[argmin])
}
