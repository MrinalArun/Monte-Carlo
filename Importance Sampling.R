# Part a - No var reduction techniques
Lmt <- 5
x <- rexp(10000,1)
y1 <- x >= Lmt

# Part b - Importance sampling
x  <- rexp(10000,1/Lmt)
y2 <- (x>=Lmt)*exp(-x)/((1/Lmt)*exp(-x/Lmt))

mean(y1)
mean(y2)
var(y1)
var(y2)