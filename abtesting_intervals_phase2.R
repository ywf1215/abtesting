# This R script implements the requirements for bug 27038 as Phase 2
# 
# Author: Kala Yang
###############################################################################


# n is the total number of unique visitors sent to the version during the test.
# values are the revenue or profit per converted unique visitor in the verion
# conf.level is the desired confidence level.
# use.approx determins whether to calculate the integral over the Mellin transform directly(when FALSE, the default), or over an approximated function (when TRUE).

# the return value is the two endpoint of the interval
auvv_ci <- function(n, values, conf.level=0.95, use.approx=FALSE) {
	stopifnot(length(values) > 3, n >= length(values), conf.level < 1, conf.level > 0.75)
	vn <- length(values)
	alpha <- 1 + vn
	beta <- 1 + (n - vn)
	m <- mean(values)
	s2 <- (sum((values - m) ^ 2)) / (vn - 1)
	dt_s <- function(x, df, shift, scale) {
		stopifnot(scale > 0)
		return(dt(x=((x - shift) / scale), df=df) * (1 / scale))
	}
	l <- vn / n
	subdivisions <- 1e+08
	rel.tol = 1e-12
	mel <- function(z) {
		stopifnot(length(z) == 1)
		tryCatch({
			c1 <- integrate(Vectorize(function(x) (1 / x) * dbeta(x, shape1=alpha, shape2=beta) * dt_s(x=(z / x), df=(vn - 1), shift=m, scale=(s2 / vn))),
					lower=1e-256, upper=1, subdivisions=subdivisions, rel.tol=rel.tol)
			c2 <- integrate(Vectorize(function(x) (1 / x) * dbeta(x, shape1=alpha, shape2=beta) * dt_s(x=(z / x), df=(vn - 1), shift=m, scale=(s2 / vn))),
					lower=l, upper=1, subdivisions=subdivisions, rel.tol=rel.tol)
			return(c1$value + c2$value)
		}, error=function(e){
			return(NA)
		})
	}
	mem <- list(zmax=NA, zmin=NA, vmax=NA, vmin=NA, init=FALSE)
	approxMel <- function(z) {
		stopifnot(length(z) == 1)
		v <- mel(z)
		if(is.na(v)) {
			if (mem$init && z > mem$zmin && z < mem$zmax) {
				v <- mem$vmin + ((mem$vmax - mem$vmin) * ((z - mem$zmin) / (mem$zmax - mem$zmin)))
			}
			if (is.na(v)) {
				step = 0.005
				z.max <- z + step
				v.max <- mel(z.max)
				while (is.na(v.max)) {
					z.max <- z.max + step
					v.max <- mel(z.max)
				}
				z.min <- z - step
				v.min <- mel(z.min)
				while (is.na(v.min)) {
					z.min <- z.min - step
					v.min <- mel(z.min)
				}
				mem <- list(zmin=z.min, zmax=z.max, vmin=v.min, vmaxj=v.max, init=TRUE)
				v <- v.min + ((v.max - v.min) * ((z - z.min) / (z.max - z.min)))
			}
		}
		stopifnot(!is.na(v))
		return(v)
	}
	e <- m * (vn / n)
	approx.zero = 1e-15
	ap.min <- e
	counter <- 1
	repeat {
		v <- approxMel(ap.min)
		if (v > approx.zero) {
			ap.min <- e - (counter * abs(e))
		} else {
			break
		}
		counter <- counter * 2
	}
	s.min <- ap.min
	s.max <- e
	repeat {
		if (s.max - s.min < 1) {
			ap.min <- s.min
			break
		}
		v <- approxMel((s.min + s.max) / 2)
		if (v > approx.zero) {
			s.max <- (s.min + s.max) / 2
		} else {
			s.min <- (s.min + s.max) / 2
		}
	}
	ap.max = e
	counter <- 1
	repeat {
		v <- approxMel(ap.max)
		if (v > approx.zero) {
			ap.max <- e + (counter * abs(e))
		} else {
			break
		}
		counter <- counter * 2
	}
	s.min <- e
	s.max <- ap.max
	repeat {
		if (s.max - s.min < 1) {
			ap.max <- s.max
			break
		}
		v <- approxMel((s.min + s.max) / 2)
		if (v > approx.zero) {
			s.min <- (s.min + s.max) / 2
		} else {
			s.max <- (s.min + s.max) / 2
		}
	}
	stopifnot(ap.min > -1e+06, ap.max < 1e+06)
	approxFunc <- NULL
	ci <- function(p, approx=FALSE) {
		stopifnot(p > 0, p < 1)
		rel.tol <- 1e-06
		if (approx && is.null(approxFunc)) {
			step <- max(0.01, min(1, (ap.max - ap.min) / 5000))
			x <- seq(from=ap.min - step, to=ap.max + step, by=step)
			y <- sapply(x, approxMel)
			approxFunc <<- approxfun(x=c(-1e+06, x, 1e+06), y=c(0, y, 0), method="linear", yleft=0, yright=0)
		}
		s.min <- NA
		s.max <- NA
		s.mid <- NA
		v <- NA
		if (p < 0.5) {
			s.min <- ap.min
			s.max <- e
		} else {
			s.min <- e
			s.max <- ap.max
		}
		repeat {
			s.mid <- (s.min + s.max) / 2
			if (p < 0.5) {
				if (approx) {
					v <- integrate(approxFunc, lower=ap.min - 1, upper=s.mid, subdivisions=subdivisions, rel.tol=rel.tol)$value
				} else {
					v <- integrate(Vectorize(approxMel), lower=ap.min - 0.01, upper=s.mid, subdivisions=subdivisions, rel.tol=rel.tol)$value
				}
			} else {
				if (approx) {
					v <- 1 - integrate(approxFunc, lower=s.mid, upper=ap.max + 1, subdivisions=subdivisions, rel.tol=rel.tol)$value
				} else {
					v <- 1 - integrate(Vectorize(approxMel), lower=s.mid, upper=ap.max + 0.01, subdivisions=subdivisions, rel.tol=rel.tol)$value
				}
			}
			if (abs(v - p) < 1e-06 || (s.max - s.min) < 0.01) {
				break
			} else if (v > p) {
				s.max <- s.mid
			} else {
				s.min <- s.mid
			}
		}
		return(list(x=s.mid, value=v, p=p))
	}
	lower <- ci((1 - conf.level) / 2, approx=use.approx)
	upper <- ci((1 + conf.level) / 2, approx=use.approx)
	return(c(lower$x, upper$x))
}

# the parameters are the endpoints for the intervals of versions A and B respectively, and the return value is the probability of version B outperforming version A
cmp <- function(a.min, a.max, b.min, b.max) {
	stopifnot(a.min <= a.max, b.min <= b.max)
	if (a.min > b.max) {
		return(0)
	}
	if (b.min > a.max) {
		return(1)
	}
	s <- (a.max - a.min) * (b.max - b.min)
	i <- (max(0, (min(a.max, b.max) - max(a.min, b.min))) ^ 2) / 2
	x <- max(0, b.min - a.min) * (b.max - b.min)
	y <- max(0, b.max - a.max) * max(0, min(a.max, b.max) - max(a.min, b.min))
	return ((i + x + y) / s)
}


# t is the time period the test had been running thus far. This can be expressed at any desired granularity -- weeks, days, hours...
# a.n and b.n are the number of unique visitors sent to versions A and B respectively over that time period.
# a.values and b.values are vectors of revenue or profit per unique visitor over the time period, for versions A and B respectively.
# threshold is the desired probability that one version would outperform the other.
# conf.level is the desired confidence level.

# the return value is an estimate of the number of additional time units left until the threshold probability would be achieved, expressed in the same time units as t
time_left <- function(t, a.n, a.values, b.n, b.values, threshold, conf.level=0.95) {
	stopifnot(t > 0, threshold > 0, threshold < 1)
	min <- 1
	max <- 1
	p <- Inf
	repeat {
		a.ci <- auvv_ci(n=(a.n * max), values=rep(x=a.values, each=max), conf.level=conf.level, use.approx=FALSE)
		b.ci <- auvv_ci(n=(b.n * max), values=rep(x=b.values, each=max), conf.level=conf.level, use.approx=FALSE)
		p <- cmp(a.min=a.ci[1], a.max=a.ci[2], b.min=b.ci[1], b.max=b.ci[2])
		if (p > threshold || p < (1 - threshold)) {
			break
		}
		min <- max
		max <- max * 2
	}
	if (max == 1) {
		return(0)
	}
	mid <- max
	repeat {
		mid <- ceiling((max + min) / 2)
		if (mid == max) {
			break
		}
		a.ci <- auvv_ci(n=(a.n * mid), values=rep(x=a.values, each=mid), conf.level=conf.level, use.approx=FALSE)
		b.ci <- auvv_ci(n=(b.n * mid), values=rep(x=b.values, each=mid), conf.level=conf.level, use.approx=FALSE)
		p <- cmp(a.min=a.ci[1], a.max=a.ci[2], b.min=b.ci[1], b.max=b.ci[2])
		if (p > threshold || p < (1 - threshold)) {
			max <- mid
		} else {
			min <- mid
		}
	}
	return((mid - 1) * t)
}

