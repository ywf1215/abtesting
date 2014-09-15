compute_pvalue <- function(Xa_hat, Xb_hat, sigma_diff, df)
{
	t <- (Xa_hat-Xb_hat)/sigma_diff
	p_value <- 2*(1-pt(abs(t),df))
	return(p_value)
}
compute_confi_intv <- function(Xa_hat, Xb_hat, sigma_diff, df, conf)
{
	alpha <- 1 - conf		# alpha for 95% confidence is 0.05.
	t_a <- qt(1-alpha/2,df)
	intv_a <- (Xa_hat-Xb_hat) - t_a*sigma_diff
	intv_b <- (Xa_hat-Xb_hat) + t_a*sigma_diff
	if (intv_a <= intv_b) {
		return(c(intv_a, intv_b))
	} else {
		return(c(inv_b, intv_a))
	}
}

# compute the p-value and confidence interval
compute_pvalue_intv <- function(Ra, Na, Rb, Nb, Va, Vb, conf=0.95, balance=FALSE)
{
	if (balance) {
		if (Na > Nb) {
			Na <- Nb
		} else {
			Nb <- Na
		}
	}
	Xa_hat <- Ra/Na
	Xb_hat <- Rb/Nb
	Sa <- sqrt(Va)
	Sb <- sqrt(Vb)
	sigma_diff <- sqrt(Sa^2/Na + Sb^2/Nb)
	df <- (Sa^2/Na+Sb^2/Nb)^2/((Sa^2/Na)^2/(Na-1)+(Sb^2/Nb)^2/(Nb-1))
	
	p_value <- compute_pvalue(Xa_hat, Xb_hat, sigma_diff, df)
	intv <- compute_confi_intv(Xa_hat, Xb_hat, sigma_diff, df, conf)
	
	return(c(p_value, intv))
}

# compute the msr(minimum sample size) for B with conf 0.95, if you want to get A's, need msr*k
compute_msr <- function(Ra, Na, Rb, Nb, Va, Vb, conf=0.95, power=0.8)
{
	Xa_hat <- Ra/Na
	Xb_hat <- Rb/Nb
	Sa <- sqrt(Va)
	Sb <- sqrt(Vb)
	k <- Na/Nb
	mu_diff <- Xa_hat - Xb_hat
	Z_alpha <- qt((1-conf)/2 + conf, Inf)			# 1.96 for 95% confidence.
	Z_beta <- qt(power, Inf)					# Z_beta is the one sided critical value for power
	msr <- ((Z_alpha+Z_beta)^2*(Sa^2/k+Sb^2))/(mu_diff^2)	# minimum sample size
	
	return(msr)
}

# compute_conf(Ra, Na, Rb, Nb, Va, Vb, msr, power=0.8)
# we pass the current volume or whatever volume for B to compute B's confidence
compute_conf <- function(Ra, Na, Rb, Nb, Va, Vb, msr, power=0.8)
{
	Xa_hat <- Ra/Na
	Xb_hat <- Rb/Nb
	Sa <- sqrt(Va)
	Sb <- sqrt(Vb)
	k <- Na/Nb
	mu_diff <- abs(Xa_hat - Xb_hat)
	rhs <-  quote({
		pnorm(mu_diff*sqrt(msr/(Sa^2/k+Sb^2))-qnorm(1-alpha/2))
	})
	
	alpha <- uniroot(function(alpha) eval(rhs) - power, c(0,2))$root

	est_conf = 1-alpha
	if(est_conf < 0){
		return(0)
	}
	return(est_conf)
}
