ker.epa = function(x,s)
{
  if(-s<x & x<s){
    return(3/4*(1-x^2/s^2))
  }else{
    return(0)
  }
}


#Stat code not correct
stat = function(n,Dz,hz,dz,Dyz,hyz,dyz,L){
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker.epa(x,hyz)),n,n)
  Fz = apply(Kz,1,sum)*(1/(n*(hz^dz)))
  Fyz = apply(Kyz,1,sum)*(1/(n*(hyz^dyz)))
  T = lapply(1:n, function(p){
    W1 = Kz[p,]/sum(Kz[p,])
    W2 = Kyz[p,]/sum(Kyz[p,])
    S1 = W1%*%t(W1)
    S2 = W2%*%t(W2)
    Q = lapply(1:n, function(i){
      S = (t(W1-W2)%*%L[[i]])^2
      return(sum(S*S1[i,]) + sum(S*S2[i,]))
    })
    return(sum(unlist(Q)))
  })
  G = unlist(T)
  G = G*Fz^4*Fyz^4
  return(sum(G)/n)
}

cBD.test = function(X, Y, Z, R=500){
  n = nrow(X)
  dz = ncol(Z)
  dy = ncol(Y)
  dyz = dy+dz

  Dx = as.matrix(dist(X))
  Dz = (dist(Z))
  Dyz = (dist(cbind(Y,Z)))
  hz = quantile(Dz,0.5)*n^(-1/(dz+2))
  hyz = quantile(Dyz,0.5)*n^(-1/(dyz+2))
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker.epa(x,hyz)),n,n)

  L = lapply(1:n, function(i){
    sapply(1:n, function(j){
      return(as.numeric(Dx[,i]<Dx[j,i]))
    })
  })
  D = stat(n,as.matrix(Dz),hz,dz,as.matrix(Dyz),hyz,dyz,L)

  #hz1 = bw.selection(X, Z, p1 =0, p2 = 0.5)

  hz1 = median(Dz)*n^(-1/3)
  Kz1 = matrix(sapply(as.matrix(Dz),function(x) ker.epa(x,hz1)),n,n)

  Pi = replicate(R,{sapply(1:n, function(i){
    W = Kz1[i,]
    if(sum(W!=0)==1){
      u = as.matrix(Dz)[i,]
      u = u[-i]
      m = min(u)
      return(which(u==m))
    }else{
      return(sample((1:n)[-i], 1, prob = W[-i]/sum(W[-i])))
    }
    }
  )}
  )
  D1 = numeric(R)
  L1 = list()

  for(i in 1:R){
    s = Pi[,i]
    for(k in 1:n){
      L1[[k]] = L[[s[k]]][s,s]
    }
    D1[i] = stat(n,as.matrix(Dz),hz,dz,as.matrix(Dyz),hyz,dyz,L1)
  }
  pval = (sum(D1==D) + sum(D1>D)+1)/(R+1)
  return(pval)
}

