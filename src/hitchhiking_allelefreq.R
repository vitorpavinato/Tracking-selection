#########################################################
# c calculation based on Durrett and Schweinsberg, 2004 #
#########################################################

cSM <- function(q0, r, s){(1-q0^(r/s))}

################################################################
# change on allele frequence xAB based on Smith and Haigh 1974 #
################################################################

xAB <- function(x0, q0, r, s){
        c = (1-q0)^(r/s)
        x = 1 - c + c*x0
        res = c(cSM = c, xABc = x, x0i = x0)
        return(res)
}

X0f <- function(xAB, q0, r, s){
        c = (1-q0)^(r/s)
        x0 = (xAB - 1 + c)/c
        res = c(cSM = c, x0c = x0, xABi = xAB)
        return(res)
}

xAB_c <- xAB(x0=0.15, q0=0.05, r=1e-7, s=0.01)
x0f_c <- X0f(xAB=0.1500004, q0=0.05, r=1e-7, s=0.01)

xAB_c[2] - x0f_c[2]
