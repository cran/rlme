dispvar <-
function (x, score = 1) 
{
    n = length(x)
    rx = rank(x, ties.method = c("random"))
    if (score == 1) {
        sc = sqrt(12) * ((rx/(n + 1)) - 0.5)
        dispvar = sqrt(pi/3) * sum(x * sc)/n
    }
    if (score == 2) {
        sc = sign((rx - (n + 1)/2))
        dispvar = sqrt(pi/2) * sum(x * sc)/n
    }
    dispvar
}
