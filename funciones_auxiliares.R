library(ks)
library(pROC)
library(here)
library(jsonlite)

GetBias = function(estimaciones , valor_real_estimador){
  bias = mean(estimaciones) - valor_real_estimador
  return(bias)
}

GetRMSE = function(estimaciones , valor_real_estimador){
  
  RMSE = sqrt(mean((estimaciones - valor_real_estimador)^2))
  return(RMSE)
}

estandarizacion = function(x) {
  return(log(x +1)/ (1 + log(x+1)) )
}

ObtenerMux = function(AUC , sigmax , sigmay , muy){
  mux = sqrt(sigmay^2 + sigmax^2) * qnorm(AUC) + muy
  return(mux)
}
ObtenerRateGamma = function(Y_pob, rate_x, auc_target, tol = 0.01, max_iter = 100){
  lower = 0
  upper = 15
  iter = 0
  
  while(iter < max_iter){
    bisectriz = (lower + upper) / 2
    iter = iter + 1
    set.seed(1)
    # Muestra de X
    X_pob = rgamma(1e5, shape = bisectriz, rate = rate_x)
    auc_iter = as.numeric(pROC::auc(response = c(rep(1, 1e5), rep(0, 1e5)), predictor = c(Y_pob, X_pob)))
    cat('AUC iteración ', iter, ': ', auc_iter, '\n')
    
    if(abs(auc_iter - auc_target) < tol){
      return(bisectriz)
    }
    if(auc_iter > auc_target){
      upper = bisectriz
    }
    if(auc_iter < auc_target){
      lower = bisectriz
    }
    
  }
  return(F)
}


# Calcula el valor de eta usando una curva ROC y su derivada en formato vector
# y una malla de integracion
eta_log_stable = function(roc, roc_prima, mesh){
  etas = numeric()
  etas[1] = ifelse(roc_prima[1] <= 1,
                   ((roc_prima[1] - 1) ^ 2) * mesh[1],
                   ((roc_prima[1] - 1) ^ 2) / roc_prima[1] * roc[1]
  )
  for (i in 2:length(mesh)) {
    etas[i] = ifelse(roc_prima[i] <= 1,
                     ((roc_prima[i] - 1) ^ 2) * (mesh[i] - mesh[i-1]),
                     ((roc_prima[i] - 1) ^ 2) / roc_prima[i] * (roc[i] - roc[i-1])
    )
  }
  eta = sum(etas)
  return(estandarizacion(eta))
}

HacerBoxCox = function(xo, yo, P =F){
  number_negative_obs = sum(xo <= 0) + 
    sum(yo <= 0)
  
  if (number_negative_obs > 0 ) {
    constant = -min(min(xo) , min(yo)) + 5e-4
    xo = xo + constant
    yo = yo + constant
  }
  likbox=function(h,data,n){
    m=length(data)-n
    x=data[1:n]
    y=data[(n+1):length(data)]
    if (abs(h)<1e-5){
      xh=log(x)
      yh=log(y)
    } else {
      xh=((x^h)-1)/h
      yh=((y^h)-1)/h
    }
    oout=-n/2*log(sum((xh-sum(xh)/n)^2)/n)-
      m/2*log(sum((yh-sum(yh)/m)^2)/m) +(h-1)*
      (sum(log(x))+sum(log(y)))
    return(-oout)
  }
  
  h_ini=-0.6
  hhat=optim(h_ini,likbox,data=c(xo,yo),n=length(xo),method="BFGS")$par
  if (abs(hhat)<1e-5){
    muestra1=log(xo)
    muestra2=log(yo)
  } else {
    muestra1=((xo^hhat)-1)/hhat
    muestra2=((yo^hhat)-1)/hhat
  }
  
  if (P) {
    cat('El lambda de la transformada es: ', hhat, '\n')
  }
  ret = list(muestra1 = muestra1, muestra2 = muestra2)
  return(ret)
}

EtaBinormal=function(controles,casos, bc = F, p=seq(0.00001,0.99999,length.out=10000)){
  
  if (bc) {
    df_trans = HacerBoxCox(controles, casos)
    controles = df_trans$muestra1
    casos = df_trans$muestra2
  }
  
  mux=mean(controles)
  muy=mean(casos)
  sigmax=sd(controles)
  sigmay=sd(casos)
  ro=sigmax/sigmay
  delta=(muy-mux)/sigmay
  ROC=1-pnorm(qnorm(1-p,mux,sigmax),mux+delta*sigmax/ro,sigmax/ro)
  ROCprima=(ro*exp(-0.5*(delta+ro*qnorm(p))^2))/exp(-0.5*(qnorm(p))^2)
  return(eta_log_stable(ROC, ROCprima, p))
}

# Obtiene la estimación de la densidad por Kernel
DensidadKernel=function(datos,puntos,h){
  ndatos=length(datos)
  npuntos=length(puntos)
  matk=dnorm((puntos%*%t(rep(1,ndatos))-t(datos%*%t(rep(1,npuntos))))/h)
  as.vector((matk%*%rep(1,ndatos))/(ndatos*h))
}

# Obtiene la estimación de la distribución por Kernel
DistKernel=function(datos,puntos,h){
  ndatos=length(datos)
  npuntos=length(puntos)
  matk=pnorm((puntos%*%t(rep(1,ndatos))-t(datos%*%t(rep(1,npuntos))))/h)
  as.vector((matk%*%rep(1,ndatos))/(ndatos))
}

# Evalua sobre la estimacion Kernel
Evaluate = function(punto, funcion, mesh){
  f_sorted = sort(funcion)
  l = length(mesh)
  posicion = sum((mesh < punto))
  r_index = 1
  if (posicion == l) {
    r_index = l -1
  }
  if ((posicion < l) & (posicion > 1)) {
    r_index = posicion
  }
  value = mean(c(funcion[r_index] ,funcion[r_index + 1]))
  return(value)
}

# Halla el cuantil de la estimacion Kernel
Inverse = function(punto, funcion, mesh){
  f_sorted = sort(funcion)
  l = length(mesh)
  posicion = sum(funcion < punto)
  r_index = 1
  if (posicion == l) {
    r_index = l -1
  }
  if ((posicion < l) & (posicion > 1)) {
    r_index = posicion
  }
  value = mean(c(mesh[r_index] ,mesh[r_index + 1]))
  return(value)
}

EtaKernel = function(muestra_sanos, muestra_enfermos, metodo, mesh_size = 1000, bc = F){
  
  if (bc) {
    df_trans = HacerBoxCox(muestra_sanos, muestra_enfermos)
    muestra_sanos = df_trans$muestra1
    muestra_enfermos = df_trans$muestra2
  }
  
  m = c(muestra_sanos, muestra_enfermos)
  if (metodo == 'optimo') {
    h = 1.06*sd(m)*length(m)^(-1/5)
    h_sanos= 1.06*sd(muestra_sanos)*length(muestra_sanos)^(-1/5)
    h_enfermos= 1.06*sd(muestra_enfermos)*length(muestra_enfermos)^(-1/5)
  }
  if (metodo == 'hscv'){
    h =  hscv(m)
    h_sanos = hscv(muestra_sanos)
    h_enfermos = hscv((muestra_enfermos)) 
  }
  
  sorted_sanos = sort(muestra_sanos)
  sorted_enfermos = sort(muestra_enfermos)
  
  mesh = seq(min(c(muestra_sanos,muestra_enfermos)),
             max(c(muestra_sanos,muestra_enfermos)),
             length.out = mesh_size)
  
  estimated_dist_sanos = DistKernel(sorted_sanos,mesh,h)
  estimated_dist_enfermos = DistKernel(sorted_enfermos,mesh,h)
  estimated_dens_sanos = DensidadKernel(sorted_sanos, mesh,h)
  estimated_dens_enfermos = DensidadKernel(sorted_enfermos, mesh,h)
  
  p = seq(0.0001,0.999, length.out = mesh_size)
  p_opp = 1-p
  
  roc = numeric(mesh_size)
  roc_prima = numeric(mesh_size)
  
  for (i in (1:mesh_size)) {
    point = p_opp[i]
    
    inv = Inverse(point, estimated_dist_sanos, mesh)
    numerador = Evaluate(inv, estimated_dens_enfermos, mesh)
    denominador = Evaluate(inv, estimated_dens_sanos, mesh)
    
    roc[i] = 1 - Evaluate(inv, estimated_dist_enfermos, mesh)
    roc_prima[i] = numerador / denominador
    
  }
  
  return(eta_log_stable(roc, roc_prima, p))
  
}

parametric_auc = function(controles, casos, bc = FALSE) {
  if (bc) {
    df_trans = HacerBoxCox(controles, casos)
    controles = df_trans$muestra1
    casos = df_trans$muestra2
  }
  mean_controles = mean(controles)
  std_controles = sd(controles)
  mean_casos = mean(casos)
  std_casos = sd(casos)
  z = (mean_casos - mean_controles) / sqrt(std_controles^2 + std_casos^2)
  auc = pnorm(z)
  return(auc)
}

parametric_youden = function(controles, casos, bc = FALSE) {
  
  if (bc) {
    df_trans = HacerBoxCox(controles, casos)
    controles = df_trans$muestra1
    casos = df_trans$muestra2
  }
  mu0 = mean(controles)
  sigma0 = sd(controles)
  mu1 = mean(casos)
  sigma1 = sd(casos)
  a = (1 / sigma0^2) - (1 / sigma1^2)
  b = 2 * ((mu1 / sigma1^2) - (mu0 / sigma0^2))
  c = (mu0^2 / sigma0^2) - (mu1^2 / sigma1^2) - 2 * log(sigma1 / sigma0)
  discriminant = b^2 - 4 * a * c
  if (discriminant < 0) {
    stop("No hay solución real para el umbral óptimo c*.")
  }
  c1 = (-b + sqrt(discriminant)) / (2 * a)
  c2 = (-b - sqrt(discriminant)) / (2 * a)
  c_opt = ifelse(c1 > min(controles) & c1 < max(casos), c1, c2)
  Se = 1 - pnorm(c_opt, mean = mu1, sd = sigma1)
  Sp = pnorm(c_opt, mean = mu0, sd = sigma0)
  J = Se + Sp - 1
  if (J > 1) {
    J_ret = 1
  } else if (J < 0) {
    J_ret = 0
  } else {
    J_ret = J
  }
  return(J_ret)
}

parametric_loc = function(controles, casos, bc = FALSE, p = seq(0.00001, 0.99999, length.out = 10000)) {
  if (bc) {
    df_trans = HacerBoxCox(controles, casos)
    controles = df_trans$muestra1
    casos = df_trans$muestra2
  }
  mux = mean(controles)
  muy = mean(casos)
  sigmax = sd(controles)
  sigmay = sd(casos)
  ro = sigmax / sigmay
  delta = (muy - mux) / sigmay
  quantiles_normal = qnorm(p)
  exp_neg_half_quantiles2 = exp(-0.5 * (quantiles_normal)^2)
  ROC = 1 - pnorm(qnorm(1 - p, mean = mux, sd = sigmax),
                  mean = mux + delta * sigmax / ro, sd = sigmax / ro)
  
  ROCprima = (ro * exp(-0.5 * (delta + ro * quantiles_normal)^2)) / 
    exp_neg_half_quantiles2
  mesh_size = length(p)
  if (ROCprima[1] <= 1) {
    loc_1 = sqrt(ROCprima[1]^2 + 1) * p[1]
  } else {
    loc_1 = sqrt(ROCprima[1]^2 + 1) / ROCprima[1] * ROC[1]
  }
  loc_rest = sapply(2:mesh_size, function(k) {
    if (ROCprima[k] <= 1) {
      sqrt(ROCprima[k]^2 + 1) * (p[k] - p[k - 1])
    } else {
      sqrt(ROCprima[k]^2 + 1) / ROCprima[k] * (ROC[k] - ROC[k - 1])
    }
  })
  locaux = c(loc_1, loc_rest)
  loc = (sum(locaux) - sqrt(2)) / (2 - sqrt(2))
  return(loc)
}
