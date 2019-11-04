# install.packages("devtools")
devtools::install_github("weltybiostat/sgpv")

library(sgpv)

##########################
###### BOMAMI Trial

## OR
sgpvalue(est.lo=0.967, est.hi=7.286, null.lo=0.9, null.hi=1.11)

## RR
sgpvalue(est.lo=0.953, est.hi=4.589, null.lo=0.9, null.hi=1.11)

## RD
sgpvalue(est.lo=0.003, est.hi=0.352, null.lo=-0.1, null.hi=0.1)
sgpvalue(est.lo=0.003, est.hi=0.352, null.lo=-0.05, null.hi=0.05)

## Logistic REgression

##########################
###### BOMAMI Trial

## Primary OR
sgpvalue(est.lo=1.18, est.hi=20.19, null.lo=0.9, null.hi=1.11)
# sgpvalue(est.lo=1.18, est.hi=20.19, null.lo=0.8, null.hi=1.2)

## Covariates
sgpvalue(est.lo=1.09, est.hi=19.28, null.lo=0.9, null.hi=1.11)
sgpvalue(est.lo=0.86, est.hi=38.29, null.lo=0.9, null.hi=1.11)

sgpvalue(est.lo=0.66, est.hi=8.98, null.lo=0.9, null.hi=1.11)
sgpvalue(est.lo=0.36, est.hi=26.92, null.lo=0.9, null.hi=1.11)

sgpvalue(est.lo=0.95, est.hi=1.11, null.lo=0.9, null.hi=1.11)
sgpvalue(est.lo=0.09, est.hi=1.59, null.lo=0.9, null.hi=1.11)

##########################
###### Leukemia Example

data(leukstats)
plotsgpv(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
		null.lo=-0.3, null.hi=0.3,
		set.order=order(leukstats$p.value),
		x.show=7000,
		plot.axis=c("TRUE","FALSE"),
		null.pt=0, outline.zone=TRUE,
		title.lab="Leukemia Example", y.lab="Fold Change (base 10)",
		x.lab="Classical p-value ranking",
		legend.on=TRUE)
axis(side=2,at=round(log(c(1/1000,1/100,1/10,1/2,1,2,10,100,1000),
	base=10),2),labels=c("1/1000","1/100","1/10","1/2",1,2,10,100,1000),
	las=2)

##########################
###### FDR rates

fdrisk(sgpval = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  std.err = 0.8,  null.weights = 'Uniform',  null.space = c(log(1/1.1), log(1.1)),  alt.weights = 'Uniform',  alt.space = 2 + c(-1,1)*qnorm(1-0.05/2)*0.8,  interval.type = 'confidence',  interval.level = 0.05)
fdrisk(sgpval = 0,  null.lo = log(1/1.1), null.hi = log(1.1),  std.err = 0.8,  null.weights = 'Point',  null.space = 0,  alt.weights = 'TruncNormal',  alt.space = 2 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
fdrisk(sgpval = 0,  null.lo = log(1/1.5), null.hi = log(1.5),  std.err = 0.8,  null.weights = 'Point',  null.space = 0,  alt.weights = 'Uniform',  alt.space = 2.5 + c(-1,1)*qnorm(1-0.041/2)*0.8,  interval.type = 'likelihood',  interval.level = 1/8)
fdrisk(sgpval = 1,  null.lo = log(1/1.5), null.hi = log(1.5),  std.err = 0.15,  null.weights = 'Uniform',  null.space = 0.01 + c(-1,1)*qnorm(1-0.041/2)*0.15,  alt.weights = 'Uniform',  alt.space = c(log(1.5), 1.25*log(1.5)),  interval.type = 'likelihood',  interval.level = 1/8)

###
##
#