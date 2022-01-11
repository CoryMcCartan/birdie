library(tidyverse)
library(tictoc)
library(rstan)
suppressMessages(library(here))

devtools::load_all(here("."))

voterfile <- here("data/nc_voters.rds")
voters = read_rds(voterfile) %>%
    mutate(id=1:n(), .before=last_name) %>%
    select(-starts_with("pred."))
#voters = make_nc_df("Duplin") # try Duplin, Hoke, Durham, Orange, Wake (largest)
p_xr = prop.table(table(voters$party, voters$race))[, c("white", "black", "hisp", "asian", "other")]
p_x = rowSums(p_xr)
p_r = colSums(p_xr)
#d = slice_sample(voters, n=1e5, replace=TRUE)
d = voters

sm = stan_model("explore/stan/model_distr.stan")

{tic()
x = model_race(party, last_name, zip, data=d, iter=500, p_r=p_r, sm=sm, verbose=TRUE)
toc()}


# marginal race probabilities
warmup = 1:100
p_r
colMeans(x$base_surn)
colMeans(x$baseline)
prop.table(table(as.integer(x$gibbs[, -warmup])))
plot(colMeans(x$gibbs == 1), type='l')

# joint race-party probabilities
p_xr_base_surn = map(rownames(p_xr), ~ colMeans(x$base_surn * (d$party == .))) %>%
    do.call(rbind, .) %>%
    `rownames<-`(rownames(p_xr))
p_xr_base = map(rownames(p_xr), ~ colMeans(x$baseline * (d$party == .))) %>%
    do.call(rbind, .) %>%
    `rownames<-`(rownames(p_xr))
p_xr_gibbs = map(1:5, function(race) {
    map_dbl(rownames(p_xr), ~ mean((x$gibbs[, -warmup] == race) * (d$party == .)))
}) %>%
    do.call(cbind, .) %>%
    `rownames<-`(rownames(p_xr)) %>%
    `colnames<-`(colnames(p_xr))
# stan
p_xr_stan =  matrix(nrow=length(p_x), ncol=length(p_r))
for (j in 1:length(p_r)) {
    for (k in 1:length(p_x)) {
        p_xr_stan[k, j] = x$stan$par[str_glue("p_xr[{k},{j}]")]
    }
}
p_xr_stan = (p_xr_stan %*% diag(p_r))  %>%
    `rownames<-`(rownames(p_xr)) %>%
    `colnames<-`(colnames(p_xr))


p_xr %*% diag(1/p_r)
p_xr_base_surn %*% diag(1/p_r)
p_xr_base %*% diag(1/p_r)
p_xr_gibbs %*% diag(1/p_r)
p_xr_stan %*% diag(1/p_r)
rake(x$lsq, p_x, p_r) %*% diag(1/p_r)
rake(x$nnls, p_x, p_r) %*% diag(1/p_r)

sum(abs(p_xr - p_xr_base_surn))
sum(abs(p_xr - rake(p_xr_base_surn, p_x, p_r)))
sum(abs(p_xr - p_xr_base))
sum(abs(p_xr - p_xr_gibbs))
sum(abs(p_xr - rake(p_xr_gibbs, p_x, p_r)))
sum(abs(p_xr - x$lsq))
sum(abs(p_xr - x$nnls))
sum(abs(p_xr - rake(x$lsq, p_x, p_r)))
sum(abs(p_xr - rake(x$nnls, p_x, p_r)))
sum(abs(p_xr - p_xr_stan))
fn_plot = Vectorize(\(a) sum(abs(p_xr - (p_xr_stan*a + (1-a)*p_xr_base_surn))))
curve(fn_plot, 0, 1)

round(p_xr_gibbs - p_xr, 3)
round(p_xr_stan - p_xr, 3)
round(p_xr_stan - p_xr_base, 3)
round(p_xr_gibbs - p_xr_base, 3)


# evaluate confounding
m = glm(I(d$party == "dem") ~ I(d$race == "white") + x$baseline[, 1], family=binomial())
#m = glm(d$lic ~ I(d$race == "white") + x$baseline[, 1], family=binomial())
summary(m)


if (F) {
    tr_ratio = function(r) {
        plot(colMeans(x[d$party == "dem", ] == r), type='l', col="blue",
             ylim=c(0, 1), ylab=names(p_r)[r])
        lines(colMeans(x[d$party == "rep", ] == r), col="red")
        abline(h=mean(d$race[d$party=="dem"] == names(p_r)[r]), col="blue", lty="dashed")
        abline(h=mean(d$race[d$party=="rep"] == names(p_r)[r]), col="red", lty="dashed")
        abline(h=p_r[r], col="grey", lty="dashed")
        abline(h=mean(d$race == names(p_r)[r]), col="grey", lty="dotted")
        lines(colMeans(x == r), col="grey")
    }
}
