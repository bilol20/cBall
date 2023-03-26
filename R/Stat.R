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



cBD = function(X,Y,Z, kernel = c("normal", "epanenchnikov")){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  Dz = dist(Z)
  Dyz = dist(cbind(Y,Z))
  hz = bw.selection(X,Z)
  hyz = bw.selection(X,cbind(Y,Z))
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

  Dx = as.matrix(dist(X))
  Dz = dist(Z)
  Dyz = dist(cbind(Y,Z))
  hz =  bw.selection(X, Z, p1 = 0.1, p2 = 0.5)
  hyz = bw.selection(X, cbind(Y,Z), p1 = 0.1, p2 = 0.5)

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

  Pi = replicate(R,{sapply(1:n, function(i) sample(1:n, 1, prob = Kz[i,]/Wz[i]))})
  D1 = resample(n,Kz,Wz,Kyz,Wyz,L,Pi)
  pval = (sum(D1>D)+1)/(R+1)
  return(pval)
}
