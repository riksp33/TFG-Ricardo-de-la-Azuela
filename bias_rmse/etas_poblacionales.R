m = 100000
l = 0.000001
u = 0.999999
eta_pob_I_n = function (mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qnorm(1-p, 0, 1)
  numerador = dnorm(inv, mux, 1)
  denominador = dnorm(inv, 0, 1)
  
  roc = 1 - pnorm(inv, mux, 1)
  roc_prima = numerador / denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_II_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qnorm(1-p, 0, 1)
  numerador = dnorm(inv, mux, 1.4)
  denominador = dnorm(inv, 0, 1)
  
  roc = 1 - pnorm(inv, mux, 1.4)
  roc_prima = numerador / denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_III_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qnorm(1-p, 0, 1)
  numerador = dnorm(inv, mux, 0.3)
  denominador = dnorm(inv, 0, 1)
  
  roc = 1 - pnorm(inv, mux, 0.3)
  roc_prima = numerador / denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}
eta_pob_IV_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qlnorm(1-p, 0, 1)
  numerador = dlnorm(inv, mux, 0.5)
  denominador = dlnorm(inv, 0, 1)
  
  roc = 1 - plnorm(inv, mux, 0.5 )
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_V_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qlnorm(1-p, 0, 1)
  numerador = dlnorm(inv, mux, 3/2)
  denominador = dlnorm(inv, 0, 1)
  
  roc = 1 - plnorm(inv, mux, 3/2 )
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_VI_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qlnorm(1-p, 0, 1)
  numerador = dlnorm(inv, mux, 0.2)
  denominador = dlnorm(inv, 0, 1)
  
  roc = 1 - plnorm(inv, mux, 0.2 )
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_VII_n = function(mux, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qlnorm(1-p, 0, 1)
  numerador = dlnorm(inv, mux, 2)
  denominador = dlnorm(inv, 0, 1)
  
  roc = 1 - plnorm(inv, mux, 2 )
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}



eta_pob_VIII_n = function(shape, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qgamma(1-p, shape = 0.5, rate = 0.5)
  numerador = dgamma(inv, shape = shape, rate = 1)
  denominador = dgamma(inv, shape = 0.5, rate = 0.5)
  
  roc = 1 - pgamma(inv, shape = shape, rate = 1)
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_IX_n = function(shape, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qgamma(1-p, shape = 0.5, rate = 0.5)
  numerador = dgamma(inv, shape = shape, rate = 4)
  denominador = dgamma(inv, shape = 0.5, rate = 0.5)
  
  roc = 1 - pgamma(inv, shape = shape, rate = 4)
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}

eta_pob_X_n = function(shape, mesh_size = m){
  p = seq(l, u, length.out = mesh_size)
  
  inv = qgamma(1-p, shape = 0.5, rate = 0.5)
  numerador = dgamma(inv, shape = shape, rate = 1/8)
  denominador = dgamma(inv, shape = 0.5, rate = 0.5)
  
  roc = 1 - pgamma(inv, shape = shape, rate = 1/8)
  roc_prima = numerador/denominador
  
  return(eta_log_stable(roc, roc_prima, p))
}