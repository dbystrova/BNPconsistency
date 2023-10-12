
#test for ratio of uniform sampling 


f_v <- function(v, K,H, n_vec, sigma){
  f_vec<- c()
  for (i in 1:K){
    res= 0
    for (j in 1:n_vec[i]){
     res= res + ((v/H)^j)*C_nk(n_vec[i], j, sigma)
    }
  }
  
  
  
}