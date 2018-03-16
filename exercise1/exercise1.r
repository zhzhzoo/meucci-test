lognormal_param <- function(expectation, variance) {
    m <- expectation;
    v <- variance;
    mu <- log(m/sqrt(1 + v/(m*m)));
    sigma2 <- log(1 + v/(m*m));
    return(list(mu=mu, sigma2=sigma2));
};

E <- 3;
V <- 5;
N <- 500;
exact <- lognormal_param(E, V);
sample <- rlnorm(N, exact$mu, sqrt(exact$sigma2));
par(mfrow=c(3, 1));

# fig 1
plot(sample);
title(sprintf("N = %d samples from lognormal(%f, %f)\nmean=%f, var=%f", N, exact$mu, exact$sigma2, mean(sample), var(sample)));

# fig 2
h <- hist(sample, breaks=50, plot=FALSE);
xfit <- seq(min(h$breaks), max(h$breaks), length = 200);
exact_pdf <- dlnorm(xfit, exact$mu, sqrt(exact$sigma2));
plot(h, ylim=c(0, max(h$density, exact_pdf)), freq=FALSE, main="");
lines(xfit, exact_pdf, col="blue");
title("Histogram of samples v.s. exact pdf");
legend("topright", "exact pdf", lty=1, col="blue");

# fig 3
empirical <- list(mu=mean(log(sample)), sigma2=var(log(sample)));
empirical_pdf <- dlnorm(xfit, empirical$mu, sqrt(empirical$sigma2));
plot(xfit, exact_pdf, type="l", lty=1, col="blue", ylim=c(0, max(exact_pdf, empirical_pdf)));
lines(xfit, empirical_pdf);
title("empirical pdf v.s. exact pdf");
legend("topright", c("exact pdf", "empirical pdf"), lty=c(1, 1), col=c("blue", "black"));
