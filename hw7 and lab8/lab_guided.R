# STA402 guided lab: partial pooling and Gibbs sampling
# All code is supplied. Run Sections 1--4 in order from the handout.
# BASE R ONLY: no package installation or internet connection is needed.
# Set the working directory to the folder containing this script and data/.
# In RStudio: Session > Set Working Directory > To Source File Location.

# SECTION 1: actual measurements and county sample sizes ------------------
if (!file.exists("data/radon_basement_mn.csv")) {
  stop("Set the working directory to the unzipped lab folder, then rerun.")
}
z <- read.csv("data/radon_basement_mn.csv", stringsAsFactors=FALSE)
y <- z$log_radon
county <- sort(unique(z$county))
g <- match(z$county, county)  # county number for each home
J <- length(county); N <- length(y)
nj <- tabulate(g, nbins=J)
ybar <- as.numeric(tapply(y, factor(g, levels=seq_len(J)), mean))
focus <- c("MURRAY", "ST LOUIS")
jfocus <- match(focus, county)
stopifnot(N==766, J==85, all(nj>0), all(is.finite(y)))
dir.create("output", showWarnings=FALSE)
raw_table <- data.frame(county=focus, n=nj[jfocus], raw_average=ybar[jfocus])
print(raw_table, row.names=FALSE, digits=4)
plot_data <- function() {
  hist(y, breaks=25, col="#DCE8EF", border="white", main="Basement radon measurements",
       xlab="Natural log of concentration (measured in pCi/L)")
}
plot_data()
png("output/data_histogram.png", width=1100, height=700, res=130)
plot_data(); dev.off()

# SECTION 2: check one conditional-mean calculation ----------------------
# First complete the two-row calculation on paper.
mu_now <- 1.5; sigma2_now <- 0.49; tau2_now <- 0.16
weight_now <- nj[jfocus]*tau2_now/(sigma2_now+nj[jfocus]*tau2_now)
mean_now <- weight_now*ybar[jfocus]+(1-weight_now)*mu_now
print(data.frame(county=focus, data_weight=weight_now,
                 conditional_mean=mean_now), row.names=FALSE, digits=4)

# SECTION 3: supplied Gibbs sampler and trace plots ----------------------
# Normal distributions use a standard deviation in rnorm().
# For v ~ IG(a,b), draw G ~ Gamma(shape=a, rate=b) and set v=1/G.
gibbs <- function(init, iter=10000, burn=2000) {
  mu <- init[1]; sigma2 <- init[2]; tau2 <- init[3]
  theta <- ybar
  saved <- matrix(NA_real_, nrow=iter-burn, ncol=J+3)
  colnames(saved) <- c("mu", "sigma2", "tau2", county)
  for (s in seq_len(iter)) {
    # Update 1: all county means, conditional on current mu and variances.
    Vtheta <- 1/(nj/sigma2 + 1/tau2)
    mtheta <- Vtheta*(nj*ybar/sigma2 + mu/tau2)
    theta <- rnorm(J, mean=mtheta, sd=sqrt(Vtheta))
    # Update 2: shared mean, using the newly sampled county means.
    Vmu <- 1/(J/tau2 + 1/2.25)
    mmu <- Vmu*(sum(theta)/tau2 + 1/2.25)
    mu <- rnorm(1, mean=mmu, sd=sqrt(Vmu))
    # Update 3: within-county variance, using the latest county means.
    sigma2 <- 1/rgamma(1, shape=2+N/2,
                       rate=1+sum((y-theta[g])^2)/2)
    # Update 4: between-county variance, using latest theta AND mu.
    tau2 <- 1/rgamma(1, shape=2+J/2,
                     rate=0.25+sum((theta-mu)^2)/2)
    # One row is the whole state after one completed sweep.
    if (s>burn) saved[s-burn,] <- c(mu, sigma2, tau2, theta)
  }
  saved
}
set.seed(4021801)
initial <- list(c(0,0.3,0.04), c(1,1,0.25), c(2,2,1), c(3,4,2))
chains <- lapply(initial, gibbs)
colors <- c("#143D59", "#168A8A", "#C16B31", "#687787")
plot_diagnostics <- function() {
  oldpar <- par(mfrow=c(2,1), mar=c(4.0,4.1,2.1,1.0), mgp=c(2.3,.7,0))
  on.exit(par(oldpar))
  for (v in c("mu","tau2")) {
    matplot(sapply(chains,function(x)x[,v]), type="l", lty=1, col=colors,
            xlab="Retained sweep", ylab=v, main=paste("Four retained chains:",v))
  }
}
plot_diagnostics()
png("output/trace_plots.png", width=1100, height=1000, res=130)
plot_diagnostics(); dev.off()
plot_acf <- function() {
  oldpar <- par(mar=c(4.5,4.1,3.6,1),mgp=c(2.3,.7,0))
  on.exit(par(oldpar))
  acf(chains[[1]][,"tau2"], main="Autocorrelation: tau2, one retained chain", lag.max=30)
}
plot_acf()
png("output/autocorrelation.png", width=1100, height=700, res=130)
plot_acf(); dev.off()
# Inspect plots before interpreting pooled summaries. Ask the TA about problems.

# SECTION 4: supplied posterior table and two-county comparison -----------
draws <- do.call(rbind, chains) # pool only after checking chain agreement
Theta <- draws[,county,drop=FALSE]
summary_draw <- function(x) c(mean=mean(x), low=unname(quantile(x,.025)),
                              high=unname(quantile(x,.975)))
county_summary <- t(sapply(focus,function(cc)summary_draw(Theta[,cc])))
rownames(county_summary) <- focus
result <- data.frame(county=focus, n=nj[jfocus], raw_average=ybar[jfocus],
                     posterior_mean=county_summary[,"mean"],
                     low=county_summary[,"low"], high=county_summary[,"high"])
print(result, row.names=FALSE, digits=4)
write.csv(result,"output/county_comparison.csv",row.names=FALSE)
plot_pooling <- function() {
  oldpar <- par(mar=c(4.2,4.4,2,1))
  on.exit(par(oldpar))
  ylim <- range(c(ybar[jfocus],county_summary,mean(draws[,"mu"])))+c(-.2,.2)
  plot(1:2, county_summary[,"mean"], xlim=c(.6,2.4), ylim=ylim,
       xaxt="n", pch=16, col="#168A8A", cex=1.3,
       xlab="County", ylab="County mean log radon", main="Raw averages and partial pooling")
  axis(1, at=1:2, labels=paste0(focus," (n=",nj[jfocus],")"))
  segments(1:2,county_summary[,"low"],1:2,county_summary[,"high"],
           col="#168A8A",lwd=3)
  points((1:2)+.08,ybar[jfocus],pch=4,cex=1.4,lwd=2,col="#687787")
  abline(h=mean(draws[,"mu"]),col="#143D59",lty=2,lwd=2)
  legend("topright",c("Posterior mean and 95% interval", "Raw county average", "Posterior mean of mu"),
         pch=c(16,4,NA),lty=c(1,NA,2),col=c("#168A8A","#687787","#143D59"),
         lwd=c(2,2,2), bty="n",cex=.8)
}
plot_pooling()
png("output/county_comparison.png",width=1100,height=800,res=130)
plot_pooling();dev.off()
cat("\nFinished. The output folder contains figures and the county table.\n")
