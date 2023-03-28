kernel = function(x,h){
  if(-h <x & x<h){
    return(0.75*(1-x^2/h^2))
  }else{
    return(0)
  }
}

fun2 = function(x) length(which(x!=0))

obj.eng = function(h,X,Y){
  n = nrow(X)
  Dx = as.matrix(dist(X))
  D = as.matrix(dist(Y))
  K = matrix(sapply(D,function(x) kernel(x,h)),n,n)
  A1 = 0
  A2 = 0
  Ind = apply(K,1,fun2)
  for(j in 1:n){
    for(i in 1:n){
      if(i!=j){
        if( Ind[j] == 1){
          y = D[j,]
          m = min(y[-j])
          f = which(y==m)
          A1 = A1 + as.integer(i==f)*Dx[i,j]
        }else{
          W = K[j,]
          A1 = A1 + K[i,j]*Dx[i,j]/sum(W[-j])
        }
      }
    }
  }

  for(i in 1:n){
    for(j in 1:n){
      for(k in 1:n){
        if(k!=i & k!=j){
          if( Ind[i] == 1 & Ind[j]==1){
              y1 = D[i,]
              m1 = min(y1[-i])
              f1 = which(y1==m1)

              y2 = D[j,]
              m2 = min(y2[-j])
              f2 = which(y2==m2)

              A2 = A2 + as.integer(k==f1)*as.integer(k==f2)*Dx[i,j]
          }
          if(Ind[i] == 1 & Ind[j]!=1){
            y1 = D[i,]
            m1 = min(y1[-i])
            f1 = which(y1==m1)

            W2 = K[j,]
            A2 = A2 + as.integer(k==f1)*Dx[i,j]*K[k,j]/sum(W2[-j])
          }
          if(Ind[j] == 1 & Ind[i]!=1){
            y1 = D[j,]
            m1 = min(y1[-j])
            f1 = which(y1==m1)

            W2 = K[i,]
            A2 = A2 + as.integer(k==f1)*Dx[i,j]*(K[i,k]/sum(W2[-i]))
          }
          if(Ind[j]!= 1 & Ind[i]!=1){

            W1 = K[i,]
            W2 = K[j,]

            A2 = A2 + K[i,k]*K[j,k]*Dx[i,j]/((sum(W1[-i])*sum(W2[-j])))
          }
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

