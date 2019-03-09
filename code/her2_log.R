her2.log<-function(prev.D, her2.LT){
# Description
# This function coverts a heritability estimate in the liability-threshold (LT) scale 
# to that in the log-risk scale.
# Arguments
# prev.D: prevalence of disease of interest
# her2.LT: heritability estimate in the LT scale
  return((dnorm(qnorm(1-prev.D))/prev.D)^2*her2.LT)
}

### Example
prev.D<-0.01
her2.LT<-0.471208  
her2<-her2.log(prev.D, her2.LT)

### AUC calculation based on heritability in the log-risk scale
pnorm(sqrt(her2/2))


### sibling relative risk 
 sqrt(exp(her2))
