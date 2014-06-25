library(cthmm)
#load the model setup and example data for a competing risks model
#observed at scheduled and informative observation times.
data(DDO_data)

#run the EM on the example data
DDO_EM=EM(rates.setup=rates.setup,
init.setup=init.setup,
emission.setup=emission.setup,
the.data=DDO_data, 
num.subjects=500,
num.states=4,
num.obs.states=3, 
tol = 1e-07, 
absorb.state=c(3,4), 
maxiter = 500, 
DDO.setup = DDO.setup,
do.DDO = T)

