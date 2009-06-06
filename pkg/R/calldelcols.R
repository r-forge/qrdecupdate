




delcols <- function(X,P=1,K=1){
  m<-nrow(X)
  ncolX<-ncol(X)
  ncolX2<-ncolX-P
  LDA<-m
  
  #input A
  qr<-qr(X)
  Q<-qr.Q(qr, complete=TRUE)
  excl <- seq(P, P+K)
  X2<-X[,-excl]
#Q_B' * C
  A<-t(Q)%*%X2
#outputs
  TAU <- vector("numeric",ncolX2)
  WORK <- vector("numeric",P+1)
  INFO<-as.integer(1)
  fort <- .Fortran("delcols", as.integer(m),as.integer(ncolX2),A= as.double(A), as.integer(LDA), as.integer(K), as.integer(P),TAU=as.double(TAU), as.double(WORK),info= INFO, PACKAGE="decUpdate")
  A <- matrix(fort$A, ncol=ncolX2)
  print(A)
  R <- A
  R[lower.tri(R)] <- 0
  if(fort$info!=0)
    warning("Error message :", fort$info)
  if(any(is.na(fort$TAU)))
    warning("TAU contains NA")
  list(A=A, TAU=fort$TAU, Qb=Q, R=R, rowX=m)
}


if(FALSE){
a <- delcols(freeny.x)
a
a$A
ma <- a$A
ma2 <- ma[1:ncol(ma),]
ma2[lower.tri(ma2)] <- 0
ma2

lower.tri(a$A)

R2 <- qr.R(qr(X2))
QR <- qr(X2)
}



####################
###DELCOLSQ
###################

delcolsq <- function(delcols,P=1,K=1){
  X <- delcols
  A <- X$A
  m<-nrow(A)
  ncolX<-ncol(A)
  ncolX2<-ncolX-P
  LDA<-m

  Q <- X$Qb
  LDQ <- m

  TAU <- X$TAU
  WORK <- vector("numeric",P+1)
  INFO<-as.integer(1)
  fort <- .Fortran("delcolsq", as.integer(m),as.integer(ncolX),A= as.double(A), as.integer(LDA),Q=as.double(Q), as.integer(LDQ), as.integer(K), as.integer(P),TAU=as.double(TAU), as.double(WORK),info= INFO, PACKAGE="decUpdate")
  if(fort$info!=0)
    warning("Error message :", fort$info)
  list(Q=t(matrix(fort$Q, nrow=delcols$rowX)), TAU=fort$TAU)
}



