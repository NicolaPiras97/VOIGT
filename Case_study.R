require(voigt)

#voigt density
dvoigt <- function(x, sigma, gamma) {
  sapply(x, function(x0) {
    integrand <- function(xp) {
      gaussian <- (1 / (sigma * sqrt(2 * pi))) * exp(-xp^2 / (2 * sigma^2))
      lorentzian <- gamma / (pi * ((x0 - xp)^2 + gamma^2))
      gaussian * lorentzian
    }
    integrate(integrand, lower = -Inf, upper = Inf,
              rel.tol = .Machine$double.eps^0.25,
              subdivisions = 1000)$value
  })
}

draman<-read.table("Raman_spectrum.txt")
x<-draman$V1
y<-draman$V2
plot(x,y,type="l", main="")

#select peak (number four)
datapeak = draman[(779+3):(1100-3),]
plot(x,y,type="l", main="")
abline(v=x[779],col="red")
abline(v=x[1100],col="red")

x<-datapeak[,1]
y<-datapeak[,2]
plot(x,y,type="l", main="")

# centering peak and scale
offset<-451.1; A<- 27294.56; medianx =407.30
x= x-medianx
y = (y-offset)/A
smooth <- smooth.spline(x, y, spar = 0.5) 
ys <- predict(smooth, x)$y

# sample voigt from histogram
int<- 1:length(x)
breaks <-  x[- length(x)] + (x[-1] - x[-length(x)] ) /2
breaks = c( x[1]- (x[2]-x[1])/2 ,breaks, x[length(x)]+ (x[length(x)]-x[length(x)-1])/2)
samplesize=5000
bins <- sample(int, size=samplesize, prob=ys, replace=TRUE)
samplehist <- vector()
for (j in 1:length(bins)){
  samplehist[ j ] = runif(1, min= breaks[ bins[ j ] ], max= breaks[ bins[ j ] +1 ])
}

# estimate using Gibbs
res <- evoigt(samplehist)
# mu, sigma and gamma point estimates
orig.par<-unlist(res["posterior mean"])
orig.par
# w_G = sqrt(8 log 2) *sigma; w_L = 2* gamma 
new.par<- orig.par*c(0, sqrt(8*log(2)),  2)
names(new.par)<-c("mu","w_G", "w_L")
new.par

# plot estimated voigt versus empirical profile (original)
plot(x+ medianx, y*A+offset,type="l",lwd=2,col="black",xlab="",ylab="")
lines(x+medianx,dvoigt(x,sigma=orig.par[2],gamma=orig.par[3])*A+offset,col="red",lwd=3,lty=1)
legend("topright",legend = c("empirical profile","fitted profile"),col = c("black", "red"),lty = c(1, 1),lwd = c(2, 3),bty = "n")

# R^2 
yhat <- dvoigt(x,sigma=orig.par[2],gamma=orig.par[3])*A+offset
my<-mean(ys*A+offset)
1- sum(ys*A+offset-yhat)^2/sum((ys*A+offset-my)^2) 
