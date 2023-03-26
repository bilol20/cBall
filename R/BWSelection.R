kernel = function(x,h){
  if(-h <x & x<h){
   return(0.75*(1-x^2/h^2))
  }else{
    return(0)
  }
}

obj.eng.split = function(h,X,Y,Z){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  V = cbind(Y,Z)
  s = sample(1:n, size = as.integer(0.55*n))
  m1 = length(s)
  m2 = n-m1
  D = as.matrix(dist(V))
  K = matrix(sapply(D,function(x) kernel(x,5)),n,n)
  K1 = K[s,-s]
  W1 = apply(K1,2,sum)
  D1 = Dx[s,-s]
  A1 = 0
  A2 = 0
  for(i in 1:m1){
    for(j in 1:m2){
      if( W1[j] == 0){
        A1 = A1+D1[i,j]/m1
      }else{
        A1 = A1 + K1[i,j]*D1[i,j]/W1[j]
      }
    }
  }

  D2 = D[s,s]
  for(i in 1:m1){
    for(j in 1:m1){
      for(k in 1:m2){
        if(W1[k]==0){
          A2 = A2 + D2[i,j]/(m1^2)
        }else{
          A2 = A2 + K1[i,k]*K1[j,k]*D2[i,j]/(W1[k]^2)
        }
      }
    }
  }

  Q = matrix(0,n,n)
  for(i in 1:n){
    Q[i,] = K[i,]/sum(K[i,])-1
  }

  T = 2*A1/(m2)-A2/(m2)+1/sum((Q)^2)

  return(T)
}

bw.selection.split = function(X,Y,Z){
  V = cbind(Y,Z)
  D = dist(V)
  q = quantile(D,prob = c(0.01,0.99))
  res = optimize(obj.eng.split,interval = q, X=X,Y=Y,Z=Z)
  return(res$minimum)
}


obj.eng.cv = function(h,X,Y){
  n = nrow(X)
  Dx = as.matrix(dist(X))

  D = as.matrix(dist(Y))
  K = matrix(sapply(D,function(x) kernel(x,h)),n,n)

  A1 = 0
  A2 = 0
  for(i in 1:n){
    r = K[i,]
    W1 = sum(r[-i])
    for(j in 1:n){
      if( W1 == 0){
        A1 = A1+D[i,j]/(n-1)
      }else{
        A1 = A1 + r[j]*D[i,j]/W1
      }
    }

  }


  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:(n-2)){
        K1 = K[i,]
        K1 = K1[-c(i,j)]
        K2 = K[j,]
        K2 = K2[-c(i,j)]
        W1 = sum(K1)
        W2 = sum(K2)
        if(W1==0&W2==0){
          A2 = A2+D[i,j]/((n-2)^2)
        }else{
          if(W1==0){
            A2 = A2+D[i,j]*K2[k]/(W2*n)
          }else{
            if(W2==0){
              A2 = A2+D[i,j]*K1[k]/(W1*n)
            }else{
              A2 = A2 + K1[k]*K2[k]*D[i,j]/(W1*W2)
            }
          }
        }
      }
    }
  }

  T = 2*A1/(n)-A2/(n)

  return(T)
}

bw.selection.cv = function(X,Y){
  D = dist(X)
  q = quantile(D,prob = c(0.1,0.5))
  res = optimize(obj.eng.cv,interval = q, X=X,Y=Y)
  return(res$minimum)
}

#Initial Thoughts
kernel = function(x,h){
  if(-h <x & x<h){
    return(0.75*(1-x^2/h^2))
  }else{
    return(0)
  }
}

fun1 = function(x){
  n = ncol(x)
  t = numeric(n)
  for(i in 1:n){
    e = x[i,]
    t[i] = sum(e[-i])
  }
  return(t)
}

fun2 = function(x) length(which(x!=0))

obj.eng = function(h,X,Y){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  D = as.matrix(dist(Y))
  K = matrix(sapply(D,function(x) kernel(x,h)),n,n)
  W = fun1(K)
  A1 = 0
  A2 = 0
  Ind = apply(K,1,fun2)
  for(i in 1:n){
    for(j in 1:n){
        if( Ind[j] == 1){
          y = D[i,]
          m = min(y[-i])
          f = which(y==m)
          A1 = A1 + Dx[i,f]
        }else{
          A1 = A1 + K[i,j]*Dx[i,j]/W[j]
        }
    }
  }

  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:n){
        if(k!=i & k!=j){
          if( Ind[k] == 1){
              A2 = A2
            }
          }else{
            W = K[k,]
            W = sum(W[-c(i,j,k)])
            A2 = A2 + K[i,k]*K[j,k]*Dx[i,j]/(W^2)
          }
        }else{
          A2 = A2
        }
      }
    }
  }

  T = 2*A1/(n)-A2/(n)
  return(T)
}

bw.selection = function(X,Y,p1,p2){
  D = dist(Y)
  if(p1>0){
    q = quantile(D,prob = c(p1,p2))
    res = optimize(obj.eng,interval = q, X = X,Y= Y)
  }else{
    q = quantile(D,prob = c(p2))
    res = optimize(obj.eng,interval = c(0,q), X = X,Y= Y)
  }
  return(res$minimum)
}

