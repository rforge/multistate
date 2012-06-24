forwardback_R <- function(x,problist,delta,emission){
        out<-.Call( "forwardback", problist, emission, x, delta)
        names(out)<-c("logalpha","logbeta","LL")
        return(out)
}

