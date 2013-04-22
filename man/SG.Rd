\name{SG}
\alias{SG}
\title{
Confidence Interval -Sison and Glaz
}
\description{
Sison and Glaz (1995) have proposed a method (SG) based on the truncated Poisson random variable and associated central and factorial moments
}
\usage{
SG(x, alpha)
}
\arguments{
 \item{x}{
x refers to the cell counts of given contingency table corresponding to a categorical data}
  \item{alpha}{
a number between 0 and 1 to get the upper 100(1-??) percentage point of the chi square distribution
}
}
\value{
lower, upper limits of multinomial proportions together with product of length of k intervals as volume of simultaneous confidence intervals }
\author{
Dr M Subbiah
}
\seealso{
\code{\link{BMDE}},\code{\link{WALD}},\code{\link{GM}}
}
\examples{
y=c(44,55,43,32,67,78)
z=0.05
SG(y,z)
function (x, alpha) 
{
    t1 = proc.time()
    sgp = function(c) {
        s = sum(x)
        k = length(x)
        b = x - c
        a = x + c
        fm1 = 0
        fm2 = 0
        fm3 = 0
        fm4 = 0
        for (i in 1:k) {
            fm1[i] = x[i] * (ppois(a[i] - 1, x[i]) - ppois(b[i] - 
                2, x[i]))/(ppois(a[i], x[i]) - ppois(b[i] - 1, 
                x[i]))
            fm2[i] = x[i]^2 * (ppois(a[i] - 2, x[i]) - ppois(b[i] - 
                3, x[i]))/(ppois(a[i], x[i]) - ppois(b[i] - 1, 
                x[i]))
            fm3[i] = x[i]^3 * (ppois(a[i] - 3, x[i]) - ppois(b[i] - 
                4, x[i]))/(ppois(a[i], x[i]) - ppois(b[i] - 1, 
                x[i]))
            fm4[i] = x[i]^4 * (ppois(a[i] - 4, x[i]) - ppois(b[i] - 
                5, x[i]))/(ppois(a[i], x[i]) - ppois(b[i] - 1, 
                x[i]))
        }
        m1 = 0
        m2 = 0
        m3 = 0
        m4 = 0
        m4t = 0
        for (i in 1:k) {
            m1[i] = fm1[i]
            m2[i] = fm2[i] + fm1[i] - (fm1[i] * fm1[i])
            m3[i] = fm3[i] + fm2[i] * (3 - (3 * fm1[i])) + (fm1[i] - 
                (3 * fm1[i] * fm1[i]) + (2 * fm1[i]^3))
            m4[i] = fm4[i] + fm3[i] * (6 - (4 * fm1[i])) + fm2[i] * 
                (7 - (12 * fm1[i]) + (6 * fm1[i]^2)) + fm1[i] - 
                (4 * fm1[i]^2) + (6 * fm1[i]^3) - (3 * fm1[i]^4)
            m4t[i] = m4[i] - (3 * m2[i]^2)
        }
        s1 = sum(m1)
        s2 = sum(m2)
        s3 = sum(m3)
        s4 = sum(m4t)
        g1 = s3/(s2^(3/2))
        g2 = s4/(s2^2)
        z = (s - s1)/sqrt(s2)
        z2 = z^2
        z3 = z^3
        z4 = z^4
        z6 = z^6
        poly = 1 + g1 * (z3 - (3 * z))/6 + g2 * (z4 - (6 * z2) + 
            3)/24 + (g1^2) * (z6 - (15 * z4) + (45 * z2) - 15)/72
        f = poly * exp(-z2/2)/sqrt(2 * pi)
        pc = 0
        for (i in 1:k) {
            pc[i] = ppois(a[i], x[i]) - ppois(b[i] - 1, x[i])
        }
        pcp = prod(pc)
        pps = 1/dpois(s, s)
        rp = pps * pcp * f/sqrt(s2)
        rp
    }
    proc.time() - t1
    t = proc.time()
    y = 0
    s = sum(x)
    M1 = 1
    M2 = s
    c = M1:M2
    M = length(c)
    for (i in 1:M) {
        y[i] = round(sgp(c[i]), 4)
    }
    j = 1
    vc = 0
    while (j <= M) {
        if (y[j] < 1 - alpha && 1 - alpha < y[j + 1]) 
            vc = j
        else vc = vc
        j = j + 1
    }
    vc
    delta = ((1 - alpha) - y[vc])/(y[vc + 1] - y[vc])
    sp = x/s
    LL = round(sp - (vc/s), 4)
    UL = round(sp + (vc/s) + (2 * delta/s), 4)
    LLA = 0
    ULA = 0
    for (r in 1:length(x)) {
        if (LL[r] < 0) 
            LLA[r] = 0
        else LLA[r] = LL[r]
        if (UL[r] > 1) 
            ULA[r] = 1
        else ULA[r] = UL[r]
    }
    diA = ULA - LLA
    VOL = round(prod(diA), 8)
    cat("Lower Limit\n")
    print(LLA)
    cat("Upper Limit\n")
    print(ULA)
    cat("Volume\n")
    print(VOL)
  }
}
\keyword{ confidencelimit }
