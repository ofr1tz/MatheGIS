x1 <- 1
y1 <- 5
x2 <- -7
y2 <- 3

s1 <- 12
s2 <- 5

s3square <- (x2-x1)^2+(y2-y1)^2
print(s3square)

p_divided_by_s3 <- 0.5*(1+(s2^2/s3square)-(s1^2/s3square))
print(p_divided_by_s3)

h_divided_by_s3 <- sqrt((s2^2/s3square)-p_divided_by_s3^2)
print(h_divided_by_s3)

x3 <- x1+p_divided_by_s3*(x2-x1)+h_divided_by_s3*(-(y2-y1))
y3 <- y1+p_divided_by_s3*(y2-y1)+h_divided_by_s3*(x2-x1)
print(c(x3, y3))

checks1 <- sqrt((x3-x1)^2+(y3-y1)^2)
print(checks1)

checks2 <- sqrt((x3-x2)^2+(y3-y2)^2)
print(checks2)