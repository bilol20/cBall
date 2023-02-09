ker = function(x,s)
{
  if(-s<x & x<s){
    return(3/4*(1-x^2/s^2))
  }else{
    return(0)
  }
}

Stat = function(X,Y,Z){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  Dz = dist(Z)
  Dyz = dist(cbind(Y,Z))
  sz = median(Dz)
  syz = median(Dyz)
  Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,sz)),n,n)
  Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker(x,syz)),n,n)
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

Test = function(X,Y,Z,R=500){
n = nrow(X)
Dx = as.matrix(dist(X))
Dz = dist(Z)
Dyz = dist(cbind(Y,Z))
sz = median(Dz)
syz = median(Dyz)
Kz = matrix(sapply(as.matrix(Dz),function(x) ker(x,sz)),n,n)
Kyz = matrix(sapply(as.matrix(Dyz),function(x) ker(x,syz)),n,n)
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
D1 = replicate(R,{
  s = numeric(n)
  for(i in 1:n){
    s[i] = sample(1:n,1,prob = Kz[i,]/Wz[i])
  }
  L1 = list()
  for(i in 1:n){
    L1[[i]] = L[[s[i]]][s,s]
  }
  CppS(n,Kz,Wz,Kyz,Wyz,L1)
  })
pval = (sum(D1>D)+1)/(R+1)
return(pval)
}
