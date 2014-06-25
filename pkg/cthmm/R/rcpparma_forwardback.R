#Changes to the forward backward code
forwardback_R <- function(x,problist,delta,emission,time.dep=F){
	if(time.dep==F){
		out<-.Call( "forwardback", problist, emission, x, delta)
	}else{
		out<-.Call( "forwardbacktimedep", problist, emission, x, delta)
	}
	names(out)<-c("logalpha","logbeta","LL")
	return(out)
}
