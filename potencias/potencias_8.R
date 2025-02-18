
set.seed(1)
B = 500
MC = 1000
alpha = 0.05
n_values = c(20, 50, 100)

results = list(title = "Resultados de potencia del test para distintos tamaños muestrales. Escenario 8")

for (n in n_values) {
  p_val_auc = numeric(MC)
  p_val_youden = numeric(MC)
  p_val_loc = numeric(MC)
  p_val_param = numeric(MC)
  p_val_kernel_opt = numeric(MC)
  p_val_kernel_hscv = numeric(MC)
  
  for (i in 1:MC) {
    
    controles = rgamma(n, shape = 0.5, rate = 0.5)
    casos = rgamma(n, shape = 1.05, rate = 1)
    auc_base = parametric_auc(controles, casos, bc = T)
    youden_base = parametric_youden(controles, casos, bc = T)
    loc_base = parametric_loc(controles, casos, bc = T)
    param_base = EtaBinormal(controles, casos, bc = T)
    kernel_opt_base = EtaKernel(controles, casos, "optimo", mesh_size = 300, bc = T)
    kernel_hscv_base = EtaKernel(controles, casos, "hscv", mesh_size = 300, bc = T)
    
    auc_bootstrap = numeric(B)
    youden_bootstrap = numeric(B)
    loc_bootstrap = numeric(B)
    binormal_bootstrap = numeric(B)
    kernel_opt_bootstrap = numeric(B)
    kernel_hscv_bootstrap = numeric(B)
    
    for (b in 1:B) {
      muestra = sample(c(casos, controles), 2 * n, replace = TRUE)
      controles_b = muestra[1:n]
      casos_b = muestra[(n + 1):(2 * n)]
      
      auc_bootstrap[b] = parametric_auc(controles_b, casos_b, bc = T)
      youden_bootstrap[b] = parametric_youden(controles_b, casos_b, bc = T)
      loc_bootstrap[b] = parametric_loc(controles_b, casos_b, bc = T)
      binormal_bootstrap[b] = EtaBinormal(controles_b, casos_b, bc = T)
      kernel_opt_bootstrap[b] = EtaKernel(controles_b, casos_b, "optimo", mesh_size = 300, bc = T)
      kernel_hscv_bootstrap[b] = EtaKernel(controles_b, casos_b, "hscv", mesh_size = 300, bc = T)
    }
    
    p_val_auc[i] = mean(auc_bootstrap >= auc_base)
    p_val_youden[i] = mean(youden_bootstrap >= youden_base)
    p_val_loc[i] = mean(loc_bootstrap >= loc_base)
    p_val_param[i] = mean(binormal_bootstrap >= param_base)
    p_val_kernel_opt[i] = mean(kernel_opt_bootstrap >= kernel_opt_base)
    p_val_kernel_hscv[i] = mean(kernel_hscv_bootstrap >= kernel_hscv_base)
  }
  
  potencia_auc = mean(p_val_auc < alpha)
  potencia_youden = mean(p_val_youden < alpha)
  potencia_loc = mean(p_val_loc < alpha)
  potencia_param = mean(p_val_param < alpha)
  potencia_kernel_opt = mean(p_val_kernel_opt < alpha)
  potencia_kernel_hscv = mean(p_val_kernel_hscv < alpha)
  
  results[[paste0("N_", n)]] = list(
    auc = potencia_auc,
    youden = potencia_youden,
    LoC = potencia_loc,
    param = potencia_param,
    kernel_opt = potencia_kernel_opt,
    kernel_hscv = potencia_kernel_hscv
  )
}

json = toJSON(results, pretty = T, digits = NA)
dir_path = here('potencias', 'jsons')
full_path = file.path(dir_path,'potencias_8.json' )
write(json, file = full_path)

