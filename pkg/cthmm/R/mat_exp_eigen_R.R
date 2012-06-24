mat_exp_eigen_R<-function(eigen_decomp,t){
      out<-.Call("mat_exp_eigen",eigen_decomp,t)
     return(out)
}

mat.exp<-function (t = 1,mat,n = 20, k = 3, method = "pade")
{
#####################################################################################
#
#
#This function is from the msm library (orginally called MatrixExp) that computes the matrix expnential,
# using the eigen decomposition method if there are no duplicated eigen values (=>matrix is diagonalizable)
# or otherwise the series approximaton method
#INPUTS: t=time interval, mat=rate matrix, n=number of items in series approximation, k=Underflow correction factor
#OUTPUTS: the transition probability matrix
#Note: need to have msm loaded
#
#
#####################################################################################
    if (!is.matrix(mat) || (nrow(mat) != ncol(mat)))
	stop("\"mat\" must be a square matrix")
    nr <- nrow(mat)
    ev <- eigen(mat)
    if (length(t) > 1)
	res <- array(dim = c(dim(mat), length(t)))
    if (any(duplicated(ev$values))) {
        for (i in seq(along = t)) {
            if (method == "series") {
                matt <- mat * t[i]/2^k
                sum <- power <- diag(nr)
                for (r in 1:n) {
					power <- matt %*% power/r
					sum <- sum + power
                }
                for (s in 1:k) sum <- sum %*% sum
                resi <- sum
            }
            else if (method == "pade") {
                resi <- .C("MatrixExpPadeR", res = double(length(mat)),
						   as.double(mat), as.integer(nr), as.double(t[i]))$res
                resi <- matrix(resi, nrow = nrow(mat))
            }
            else stop("Method should be \"pade\" or \"series\"")
            if (length(t) == 1)
			res <- resi
            else res[, , i] <- resi
        }
    }
    else {
        evinv <- solve(ev$vectors)
        for (i in seq(along = t)) {
            resi <- ev$vectors %*% diag(exp(ev$values * t[i])) %*%
			evinv
            if (length(t) == 1)
			res <- resi
            else res[, , i] <- resi
        }
    }
    res
}


 