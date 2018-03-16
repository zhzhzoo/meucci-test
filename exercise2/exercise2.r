is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' maximum likelihood estimates for multivariate Student-t distribution
#'
#' @param sample    A matrix. Each column represents an observation
#' @param weight    A vector. Relative weights of the observations, its length should match number of rows of sample
#' @param nu        An integer. Degrees of freedom.
#' @param epsilon   A double. Convergence threshold.
#'
#' @return mu       Location parameter.
#' @return sigma_2  Dispersion parameter.
fit_locdisp_mlf <- function(sample, weight, nu, epsilon) {
    d = dim(e);
    if (length(d) != 2) {
        stop('The sample must be a matrix');
    }
    i = d[1];
    t = d[2];
    if (t != length(weight)) {
        stop('The number of samples and weights must match');
    }
    if (!is.wholenumber(nu) || nu < 0.5) {
        stop('The argument "nu" must be a positive integer');
    }

    # initialization
    mu = sample %*% weight;
    sigma_2 = p[1] * (sample[,1] - mu) %*% t(sample[,1] - mu);
    for (k in 2:t) {
        sigma_2 = sigma_2 + p[k] * (sample[,k] - mu) %*% t(sample[,k] - mu);
    }
    if (nu > 2) {
        sigma_2 = sigma_2 *((nu - 2) / nu);
    }

    # iteration
    prev_mu = mu;
    prev_sigma_2 = sigma_2;
    w = double(t);
    q = double(t);
    repeat {
        sigma_2_inv = solve(sigma_2);
        for (k in 1:t) {
            w[k] = (nu + i) / (nu + t(sample[,k] - mu) %*% sigma_2_inv %*% (sample[,k] - mu));
            q[k] = p[k] * w[k];
        }
        q = q / sum(q);
        prev_mu = mu;
        prev_sigma_2 = sigma_2;
        mu = sample %*% q;
        sigma_2 = p[1] * (sample[,1] - mu) %*% t(sample[,1] - mu);
        for (k in 2:t) {
            sigma_2 = sigma_2 + p[k] * (sample[,k] - mu) %*% t(sample[,k] - mu);
        }
        if (norm(prev_mu - mu, type="F") / norm(mu, type="F") < epsilon
            && norm(prev_sigma_2 - sigma_2, type="F") / norm(sigma_2, type="F") < epsilon) {
            return(list(mu=mu, sigma_2=sigma_2));
        }
    }
};

# test
t = 1000;
nu = 100;
epsilon = 1e-9;
e = matrix(rnorm(t * 2, 0, 1), 2, t);
p = rep(1 / t, t);
print(rowMeans(e));
print(cov(t(e)));
fit_locdisp_mlf(e, p, nu, epsilon);
