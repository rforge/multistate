library(cthmm)
#load the model setup and example data for a competing risks model
#observed at scheduled and informative observation times.
data(DDO_data)

covariance_DDO=get_covariance(par=DDO_EM$param,
the.data=DDO_data,
num.subjects=500,
num.states=4,
num.obs.states=3,
rates.setup=rates.setup,
emission.setup=emission.setup,
init.setup=init.setup,
DDO.setup=DDO.setup,
do.DDO=T)