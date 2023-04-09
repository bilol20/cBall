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

  s = sample(1:n, size = as.integer(0.5*n))
  Xtr = X[s,]
  Ytr = Y[s,]
  Ztr = Z[s,]

  Xte = X[-s,]
  Yte = Y[-s,]
  Zte = Z[-s,]

  n1 = length(s)
  n2 = n-n1
  hz = bw.selection(Xtr, Ztr, p1 = 0.4, p2 = 0.9)
  hyz = bw.selection(Xtr, cbind(Ytr,Ztr), p1 = 0.3, p2 = 0.7)

  Dx = as.matrix(dist(Xte))
  Dz = dist(Zte)
  Dyz = dist(cbind(Yte,Zte))
  #hz =  quantile(Dz,prob = 0.5)
  #hyz = quantile(Dyz,prob = 0.5)

  Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,hz)),n2,n2)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker(x,hyz)),n2,n2)
  Wz = apply(Kz,1,sum)
  Wyz = apply(Kyz,1,sum)
  L = lapply(1:n2, function(i){
    sapply(1:n2, function(j){
      sapply(1:n2, function(k){
        return(as.numeric(Dx[k,i]<Dx[j,i]))
      })
    })
  })
  D = CppS(n2,Kz,Wz,Kyz,Wyz,L)

  hz1 =  bw.selection(Xtr, Ztr, p1 = 0.05, p2 = 0.4)
  Kz1 = matrix(sapply(as.matrix(Dz),function(x) ker(x,hz1)),n2,n2)
  Wz1 = apply(Kz1,1,sum)


  Pi = replicate(R,{sapply(1:n2, function(i){
    W = Kz1[i,]
    sample((1:n2), 1, prob = W/sum(W))}
  )}
  )
  D1 = resample(n2,Kz,Wz,Kyz,Wyz,L,Pi)
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}
