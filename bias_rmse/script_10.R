set.seed(1)
ns = c(20, 50, 100)
AUCs = c(0.6, 0.75, 0.9)
MC = 1000
lista_json = list()
lista_json['header'] = paste('Simulacion para nueva población 10', "Tamaño Mc: ", MC)

for (auc in AUCs) {
  resultado_auc = list()
  
  estimaciones_param = numeric(MC)
  estimaciones_kernel_h_opt = numeric(MC)
  estimaciones_kernel_h_fun = numeric(MC)
  
  
  #Población sana
  shape_y = 0.5
  rate_y = 0.5
  Y_pob = rgamma(1e5, shape = shape_y, rate = rate_y)
  
  # Población enferma
  rate_x = 1/8
  shape_x = ObtenerRateGamma(Y_pob, rate_x, auc_target = auc)
  X_pob = rgamma(1e5, shape = shape_x, rate = rate_x)
  
  # Eta poblacional
  eta_true = eta_pob_X_n(shape_x)
  
  # AUC observado
  obs_auc = as.numeric(pROC::auc(response = c(rep(1, 1e5), rep(0, 1e5)),
                                 predictor = c(Y_pob, X_pob)))
  stopifnot(abs(obs_auc - auc) <= 0.05)
  
  cat('===================================================================\n')
  cat('El AUC OBSERVADO es: ', obs_auc, '\n')
  cat('El eta_true para AUC = ' , auc, 'es:', eta_true, '\n')
  cat('===================================================================\n')
  
  for (n in ns) {
    for (i in 1:MC) {
      # Poblaciones del montecarlo
      X = rgamma(n, shape = shape_x, rate = rate_x)
      Y = rgamma(n, shape = shape_y, rate = rate_y)
      
      # Vectores de estimaciones
      estimaciones_param[i] = EtaBinormal(casos = X, controles = Y, bc= T)
      estimaciones_kernel_h_opt[i] = EtaKernel(muestra_sanos= Y, muestra_enfermos= X, metodo = 'optimo', mesh_size = 300)
      estimaciones_kernel_h_fun[i] = EtaKernel(muestra_sanos= Y, muestra_enfermos= X, metodo = 'hscv', mesh_size = 300)
    }
    cat('n = ' , n, ', AUC = ', auc, 'completado')
    
    size = paste0('size:', n)
    resultado_auc[[size]] = list(
      eta_pob = eta_true,
      param = list(
        bias = GetBias(estimaciones_param, eta_true),
        rmse = GetRMSE(estimaciones_param, eta_true)
      ),
      kernel_opt = list(
        bias = GetBias(estimaciones_kernel_h_opt, eta_true),
        rmse = GetRMSE(estimaciones_kernel_h_opt, eta_true)
      ),
      kernel_hscv = list(
        bias = GetBias(estimaciones_kernel_h_fun, eta_true),
        rmse = GetRMSE(estimaciones_kernel_h_fun, eta_true)
      )
    )
  }
  lista_json[[paste0('AUC:', auc)]] = resultado_auc
}
json = toJSON(lista_json, pretty = T, digits = NA)
dir_path = here('bias_rmse', 'jsons')
full_path = file.path(dir_path,'tabla_10.json' )
write(json, file = full_path)