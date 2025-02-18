set.seed(1)
ns = c(20, 50, 100)
AUCs = c(0.6, 0.75, 0.9)
MC = 1000
lista_json = list()
lista_json['header'] = paste('Simulacion para nueva población 7', "Tamaño Mc: ", MC)

for (auc in AUCs) {
  resultado_auc = list()
  
  estimaciones_param = numeric(MC)
  estimaciones_kernel_h_opt = numeric(MC)
  estimaciones_kernel_h_fun = numeric(MC)
  
  #Población sana
  mu_y = 0
  sd_y = 1
  Y_pob = rlnorm(1e5, mu_y, sd_y)
  
  # Población enferma
  sd_x = 2
  mu_x = ObtenerMux(auc, sd_x, sd_y, mu_y)
  X_pob = rlnorm(1e5, mu_x, sd_x)
  
  # Eta poblacional
  eta_true = eta_pob_VII_n(mu_x)
  
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
      X = rlnorm(n, mu_x, sd_x)
      Y = rlnorm(n, mu_y, sd_y)
      
      # Vectores de estimaciones
      estimaciones_param[i] = EtaBinormal(casos = X, controles = Y, bc= T)
      estimaciones_kernel_h_opt[i] = EtaKernel(muestra_sanos= Y, muestra_enfermos= X, metodo = 'optimo', mesh_size = 300, bc = T)
      estimaciones_kernel_h_fun[i] = EtaKernel(muestra_sanos= Y, muestra_enfermos= X, metodo = 'hscv', mesh_size = 300, bc = T)
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
full_path = file.path(dir_path,'tabla_7.json' )
write(json, file = full_path)


