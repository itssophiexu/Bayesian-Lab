# STA402: Gaussian hierarchical radon model. Run from this folder.
# Base R; install coda once if needed: install.packages("coda")
if (!requireNamespace("coda", quietly=TRUE)) stop("Install coda: install.packages('coda')")
z <- read.csv("data/radon_basement_mn.csv", stringsAsFactors=FALSE)
y <- z$log_radon
county <- sort(unique(z$county)); g <- match(z$county, county)
J <- length(county); N <- length(y)
nj <- tabulate(g, nbins=J)
ybar <- as.numeric(tapply(y, factor(g, levels=seq_len(J)), mean))
stopifnot(N==766, J==85, all(nj>0), all(is.finite(y)))
print(range(nj)); print(c(MURRAY=nj[match("MURRAY",county)], ST_LOUIS=nj[match("ST LOUIS",county)]))
hist(y, breaks=25, xlab="Natural log radon (concentration in pCi/L)", main="Basement measurements")

gibbs <- function(init, iter=10000, burn=2000, tau_scale=0.25) {
  mu <- init[1]; sigma2 <- init[2]; tau2 <- init[3]
  theta <- ybar
  out <- matrix(NA_real_, iter-burn, J+3)
  colnames(out) <- c("mu", "sigma2", "tau2", county)
  for (s in seq_len(iter)) {
    Vtheta <- 1/(nj/sigma2 + 1/tau2)
    mtheta <- Vtheta*(nj*ybar/sigma2 + mu/tau2)
    theta <- rnorm(J, mean=mtheta, sd=sqrt(Vtheta))
    Vmu <- 1/(J/tau2 + 1/2.25)
    mmu <- Vmu*(sum(theta)/tau2 + 1/2.25)
    mu <- rnorm(1, mmu, sqrt(Vmu))
    sigma2 <- 1/rgamma(1, shape=2+N/2,
                       rate=1+sum((y-theta[g])^2)/2)

    tau2 <- 1/rgamma(1, shape=2+J/2, rate=tau_scale+sum((theta-mu)^2)/2)
    if (s>burn) out[s-burn,] <- c(mu, sigma2, tau2, theta)
  }
  out
}
set.seed(4021801)
initial <- list(c(0,0.3,0.04), c(1,1,0.25), c(2,2,1), c(3,4,2))
chains <- lapply(initial, gibbs)
ml <- do.call(coda::mcmc.list, lapply(chains, coda::mcmc))
ess <- coda::effectiveSize(ml)
# Conventional variance-based split Rhat is supplied as a supplementary check.
split <- unlist(lapply(chains, function(a) {
  h <- nrow(a)%/%2
  list(coda::mcmc(a[1:h,]), coda::mcmc(a[(h+1):(2*h),]))
}), recursive=FALSE)
rhat <- coda::gelman.diag(do.call(coda::mcmc.list,split),
                         autoburnin=FALSE, multivariate=FALSE)$psrf[,1]
print(cbind(ESS=ess[1:3], split_Rhat=rhat[1:3]))
par(mfrow=c(3,2), mar=c(3,3.5,2,1))
cols <- c("#143D59", "#168A8A", "#C16B31", "#687787")
for (v in c("mu","sigma2","tau2")) {
  matplot(sapply(chains,function(x)x[1:500,v]), type="l", lty=1,
          col=cols, xlab="Retained sweep", ylab=v)
  acf(chains[[1]][,v], main=paste("ACF:",v))
}
par(mfrow=c(1,1))
# Pool only after inspecting diagnostics.
draws <- do.call(rbind,chains); S <- nrow(draws)
Theta <- draws[,county,drop=FALSE]
mu <- draws[,"mu"]; sigma2 <- draws[,"sigma2"]; tau2 <- draws[,"tau2"]
summary_draw <- function(x) c(mean=mean(x), q025=unname(quantile(x,.025)),
                              q975=unname(quantile(x,.975)))
print(rbind(mu=summary_draw(mu),sigma=summary_draw(sqrt(sigma2)),
            tau=summary_draw(sqrt(tau2))))
# Pooling plot: add posterior county intervals to your report table separately.
plot(ybar,colMeans(Theta), pch=16, cex=.6+sqrt(nj)/8,
     xlab="Raw county mean log radon", ylab="Posterior mean county log radon")
abline(0,1,col="gray50",lty=2); abline(h=mean(mu),col="#168A8A",lwd=2)
legend("topleft",c("No-pooling reference", "Posterior mean of mu"),
       lty=c(2,1),col=c("gray50","#168A8A"),bty="n",cex=.8)
# TODO: produce the requested county and original-scale summaries.
# Access an entire county's retained theta draws with Theta[, "MURRAY"].

# TODO: for MURRAY, estimate the two distinct probabilities and simulate
# one new home per retained state. Store new_home as ORIGINAL-SCALE values.
new_home <- NEW_HOME

# Fresh measurements in the SAME represented counties, not fresh counties.
set.seed(4021802)
selected <- sample.int(S, 2000)
Yrep <- matrix(NA_real_, length(selected), N)
for (s in seq_along(selected)) {
  k <- selected[s]
  # TODO: generate N log measurements using the selected county means and sigma.
  Yrep[s,] <- REPLICATED_MEASUREMENTS
}
# TODO: one minimum per ROW, then the diagnostic plot and summaries.
