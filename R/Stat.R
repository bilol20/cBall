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



cBD = function(X,Y,Z, hz, hyz, kernel = c("normal", "epanenchnikov")){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  Dz = dist(Z)
  Dyz = dist(cbind(Y,Z))
  if(kernel == "normal"){
    ker = ker.normal
  }else{
    if(kernel == "epanenchnikov"){
      ker = ker.epa
    }else{
      print("Please specify the kernel correctly.")
    }
  }
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,hz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker(x,hyz)),n,n)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  L = lapply(1:n, function(i){
    sapply(1:n, function(j){
      sapply(1:n, function(k){
        return(as.numeric(Dx[k,i]<Dx[j,i]))
      })
    })
  })
  D = CppS(n,Kz,Wz,Kyz,Wyz,L)
  return(D)
}



cBD.test = function(X, Y, Z, R=500, beta = 0.3, kernel = c("normal", "epanenchnikov")){
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

  if(class(X)[1]=="matrix"){
    n = nrow(X)
  }else{
    n = length(X)
  }

  s = sample(n, size = as.integer(beta*n))
  n1 = length(s)
  hz =  bw.selection(X[s,], Z[s,], p1 = 0.01, p2 = 0.5)
  hyz = bw.selection(X[s,], cbind(Y,Z)[s,], p1 = 0.01, p2 = 0.5)

  Dx = as.matrix(dist(X[-s,]))
  Dz = as.matrix(dist(Z[-s,]))
  Dyz = as.matrix(dist(cbind(Y,Z)[-s,]))
  Kz = matrix(sapply(Dz,function(x) ker(x,hz)),(n-n1),(n-n1))
  Kyz = matrix(sapply(Dyz,function(x) ker(x,hyz)),(n-n1),(n-n1))
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  L = lapply(1:(n-n1), function(i){
    sapply(1:(n-n1), function(j){
      sapply(1:(n-n1), function(k){
        return(as.numeric(Dx[k,i]<Dx[j,i]))
      })
    })
  })

  D = CppS((n-n1),Kz,Wz,Kyz,Wyz,L)
  D1 = numeric(R)
  for(i in 1:R){
    Pi = sapply(1:(n-n1), function(l){
      W = Kz[l,]
      if(length(which(W!=0))==0){
        z = Dz[l,]
        m = min(z[-l])
        f = which(z==m)
        return(f)
      }else{
        return(sample(1:(n-n1), 1, prob = W/sum(W)))}
    }
    )
    L1 = list()
    for(q in 1:(n-n1)){
      L1[[q]] = L[[Pi[q]]][Pi,Pi]
    }
    D1[i] =  CppS((n-n1),Kz,Wz,Kyz,Wyz,L1)
  }
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}

