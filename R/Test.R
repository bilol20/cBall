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
stat = function(n,Dz,hz,Dyz,hyz,L, kernel = c("normal", "epanenchnikov")){
  if(kernel == "normal"){
    ker = ker.normal
  }else{
    if(kernel == "epanenchnikov"){
      ker = ker.epa
    }else{
      print("Please specify the kernel correctly.")
      exit
    }
  }
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker(x,hyz)),n,n)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  T = 0
  for(p in 1:n){
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


    for(i in (1:n)[-p]){
      for(j in (1:n)[-p]){
        T = T + sum((W1-W2)*L[[i]][,j])^2*(W1[i]*W1[j]+W2[i]*W2[j])
      }
    }
  }
  return(T/n)
}


#calibration part is not clear
cBD.test = function(X, Y, Z, beta=0.5, R=500, kernel = c("normal", "epanenchnikov")){
  if(class(X)[1]=="matrix"){
    n = nrow(X)
  }else{
    n = length(X)
  }

  s = sample(1:n, size = as.integer(beta*n))
  Xtr = X[s,]
  Ytr = Y[s,]
  Ztr = Z[s,]

  Xte = X[-s,]
  Yte = Y[-s,]
  Zte = Z[-s,]

  n1 = length(s)
  n2 = n-n1
  hz = bw.selection(Xtr, Ztr, p1 = 0.05, p2 = 0.5)
  hyz = bw.selection(Xtr, cbind(Ytr,Ztr), p1 = 0.05, p2 = 0.5)

  Dx = as.matrix(dist(Xte))
  Dz = as.matrix(dist(Zte))
  Dyz = as.matrix(dist(cbind(Yte,Zte)))
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,hz)),n2,n2)

  L = lapply(1:n2, function(i){
    sapply(1:n2, function(j){
      sapply(1:n2, function(k){
        return(as.numeric(Dx[k,i]<Dx[j,i]))
      })
    })
  })
  D = stat(n2,Dz,hz,Dyz,hyz,L, kernel)
  D1 = numeric(R)
  d = ncol(Xte)
  X.f = matrix(0, ncol = d, nrow = n2)
  re = matrix(0, ncol = d, nrow = n2)
  for(i in 1:n2){
    if(length(which(Kz[i,]!=0))>1){
      W = Kz[i,]
      for(j in (1:n2)[-i]){
        X.f[i,] = X.f[i,] + W[j]*Xte[j,]/sum(W[-i])
      }
    }else{
      x = Dz[i,]
      x1 = x[-i]
      m = min(x1)
      ff = which(x==m)
      X.f[i,] = Xte[ff,]
    }
    re[i,] = Xte[i,]-X.f[i,]
  }
  for(i in 1:R){
    s = sample(n2)
    X1 = X.f + re[s,]
    Dx1 = as.matrix(dist(X1))
    L1 = lapply(1:n2, function(i){
      sapply(1:n2, function(j){
        sapply(1:n2, function(k){
          return(as.numeric(Dx1[k,i]<Dx1[j,i]))
        })
      })
    })
    D1[i] = stat(n2,Dz,hz,Dyz,hyz,L1, kernel)
  }
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}
