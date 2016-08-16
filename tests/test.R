
x <- 1
h <- function() {
	y <- 2
	i <- function() {
		z <- 3
		c(x, y, z)
	}
	i()
}
y
h()