ker.epa = function(x,s)
{
  if(-s<x & x<s){
    return(3/4*(1-x^2/s^2))
  }else{
    return(0)
  }
}

ker.normal = function(x,s)
{
  return(dnorm(x,0,s))
}


#Stat code correct
stat = function(n,Dz,hz,Dyz,hyz,L){
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker.epa(x,hyz)),n,n)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  T = lapply(1:n, function(p){
    W1 = rep(0,n)
    W2 = rep(0,n)
    if(length(which(Kz[p,]!=0))>1&length(which(Kyz[p,]!=0))>1){
      W1 = Kz[p,]/Wz[p]
      W2 = Kyz[p,]/Wyz[p]
    }
    if(length(which(Kz[p,]!=0))==1&length(which(Kyz[p,]!=0))>1){
      x = Dz[p,]
      x1 = x[-p]
      m = min(x1)
      ff = which(x==m)
      W1[ff] = 1
      W2 = Kyz[p,]/Wyz[p]
    }
    if(length(which(Kz[p,]!=0))>1&length(which(Kyz[p,]!=0))==1){
      x = Dyz[p,]
      x1 = x[-p]
      m = min(x1)
      ff = which(x==m)
      W2[ff] = 1
      W1 = Kz[p,]/Wz[p]
    }
    if(length(which(Kz[p,]!=0))==1&length(which(Kyz[p,]!=0))==1){
      x = Dz[p,]
      x1 = x[-p]
      m = min(x1)
      ff = which(x==m)
      W1[ff] = 1
      x = Dyz[p,]
      x1 = x[-p]
      m = min(x1)
      ff = which(x==m)
      W2[ff] = 1
    }

    S1 = W1%*%t(W1)
    S2 = W2%*%t(W2)
    Q = lapply(1:n, function(i){
      S = (t(W1-W2)%*%L[[i]])^2
      return(sum(S*S1[i,]) + sum(S*S2[i,]))
    })
    return(sum(unlist(Q)))
  })
  return(sum(unlist(T)/n))
}

#calibration part is not clear
cBD.test = function(X, Y, Z, R=500){
  n = nrow(X)
  dz = ncol(Z)
  dy = ncol(Y)
  dyz = dy+dz

  Dx = as.matrix(dist(X))
  Dz = (dist(Z))
  Dyz = (dist(cbind(Y,Z)))
  hz = median(Dz)*n^(-1/(dz+2))
  hyz =median(Dyz)*n^(-1/(dyz+2))
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker.epa(x,hyz)),n,n)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)

  L = lapply(1:n, function(i){
    sapply(1:n, function(j){
      return(as.numeric(Dx[,i]<Dx[j,i]))
    })
  })
  D = stat(n,as.matrix(Dz),hz,as.matrix(Dyz),hyz,L)


  hz1 =  quantile(Dz, 0.1)* n^(-1/(dz+2))
  Kz1 = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz1)),n,n)
  Wz1 = apply(Kz1,1,sum)


  Pi = replicate(R,{sapply(1:n, function(i){
    W = Kz1[i,]
    sample((1:n), 1, prob = W/sum(W))}
  )}
  )
  D1 = numeric(R)
  for(i in 1:R){
    s = Pi[,i]
    X1 = X[s,]
    Dx1 = as.matrix(dist(X1))
    L1 = lapply(1:n, function(i){
      sapply(1:n, function(j){
        return(as.numeric(Dx1[,i]<Dx1[j,i]))
      })
    })
    D1[i] = stat(n,as.matrix(Dz),hz,as.matrix(Dyz),hyz,L1)
  }
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}

