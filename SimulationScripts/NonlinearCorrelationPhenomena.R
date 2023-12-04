# same function as below: linear interpolation between t and abs(t)
nlX <- function(al, t){(1-al)*t+al*abs(t)}

# function to induce nonlinearity in H->Y
# linear interpolation between t and abs(t)
nlY <- function(bet, t){(1-bet)*t+bet*abs(t)}


# fix n = 400, p = 500

q <- 5
n <-  400
p <- 500


Gamma <- matrix(runif(q*p, min = -1, max = 1), nrow=q)
psi <- runif(q, min = 0, max = 2)
H <- matrix(rnorm(n*q), nrow = n)
E <- matrix(rnorm(n*p), nrow = n)
e <- 0.5*rnorm(n)

al <- 1
bet <- 0

X <- nlX(al, H %*% Gamma)
Y <- nlY(bet, H %*% psi)

cor.vec <- apply(X, 2, fun <- function(v){cor(v, Y)})

plot(cor.vec)
hist(cor.vec)

plot(svd(X)$d)

