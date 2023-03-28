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

cBD.test = function(X, Y, Z, R=500, kernel = c("normal", "epanenchnikov")){
  if(class(X)[1]=="matrix"){
    n = nrow(X)
  }else{
    n = length(X)
  }
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

  s = sample(n, size = as.integer(0.5*n))
  Xtr = X[s,]
  Ytr = Y[s,]
  Ztr = Z[s,]

  Xte = X[-s,]
  Yte = Y[-s,]
  Zte = Z[-s,]

  hz =  bw.selection(Xtr, Ztr, p1 = 0.01, p2 = 0.5)
  hyz = bw.selection(Xtr, cbind(Ytr,Ztr), p1 = 0.01, p2 = 0.5)

  Dx = as.matrix(dist(Xte))
  n1 = ncol(Dx)
  Dz = as.matrix(dist(Zte))
  Dyz = as.matrix(dist(cbind(Yte,Zte)))
  Kz = matrix(sapply(Dz,function(x) ker(x,hz)),n1,n1)
  Kyz = matrix(sapply(Dyz,function(x) ker(x,hyz)),n1,n1)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  L = lapply(1:n1, function(i){
    sapply(1:n1, function(j){
      sapply(1:n1, function(k){
        return(as.numeric(Dx[k,i]<Dx[j,i]))
      })
    })
  })
  D = CppS(n1,Kz,Wz,Kyz,Wyz,L)

  Pi = replicate(R,{sapply(1:n1, function(i){
    W = Kz[i,]
    if(length(which(W!=0))==1){
      z = Dz[i,]
      m = min(z[-i])
      f = which(z==m)
      return(f)
    }else{
      sample((1:n1)[-i], 1, prob = W[-i]/sum(W[-i]))}
    }
  )}
  )

  D1 = resample(n1,Kz,Wz,Kyz,Wyz,L,Pi)
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}

