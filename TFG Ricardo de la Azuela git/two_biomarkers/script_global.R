library(MASS)

scenario_n_two_biomarkers = function(m_1_enfermos, m_2_enfermos, v_1_enfermos, v_2_enfermos, boxcox, num)
{
  set.seed(1)
  B = 500
  MC = 1000
  alpha = 0.05
  n_values = c(20, 50, 100)
  
  sigma_1_cov = sqrt(v_1_enfermos)
  sigma_2_cov = sqrt(v_2_enfermos)
  covariance_enfermos_05 = 0.5 * sigma_1_cov * sigma_2_cov
  
  covariance_enfermos_values = c(0.0, covariance_enfermos_05)
  
  results = list(title = paste("Resultados de potencia del test con múltiples marcadores para distintos tamaños muestrales. Escenario", num))
  
  for (covariance_enfermos in covariance_enfermos_values) {
    covariance_results = list()
    
    for (n in n_values) {
      p_val_auc = numeric(MC)
      p_val_youden = numeric(MC)
      p_val_loc = numeric(MC)
      p_val_param = numeric(MC)
      p_val_kernel_opt = numeric(MC)
      p_val_kernel_hscv = numeric(MC)
      
      m_1_sanos = 0.0
      m_2_sanos = 0.0
      
      v_1_sanos = 1.0
      v_2_sanos = 1.0
      
      covariance_sanos = 0.0
      
      mean_vector_sanos = c(m_1_sanos, m_2_sanos)
      mean_vector_enfermos = c(m_1_enfermos, m_2_enfermos)
      
      cov_matrix_sanos = matrix(c(v_1_sanos, covariance_sanos, 
                                  covariance_sanos, v_2_sanos), nrow = 2)
      
      cov_matrix_enfermos = matrix(c(v_1_enfermos, covariance_enfermos, 
                                     covariance_enfermos, v_2_enfermos), nrow = 2)
      
      for (i in 1:MC) {
        # Muestras de sanos y enfermos
        muestra_mv_sanos = mvrnorm(n, mean_vector_sanos, cov_matrix_sanos)
        sanos_marker_1 = muestra_mv_sanos[, 1]
        sanos_marker_2 = muestra_mv_sanos[, 2]
        
        muestra_mv_enfermos = mvrnorm(n, mean_vector_enfermos, cov_matrix_enfermos)
        enfermos_marker_1 = muestra_mv_enfermos[, 1]
        enfermos_marker_2 = muestra_mv_enfermos[, 2]
        
        # Inicialización de los vectores de bootstrap
        auc_bootstrap_marker_1 = numeric(B)
        auc_bootstrap_marker_2 = numeric(B)
        youden_bootstrap_marker_1 = numeric(B)
        youden_bootstrap_marker_2 = numeric(B)
        binormal_bootstrap_marker_1 = numeric(B)
        binormal_bootstrap_marker_2 = numeric(B)
        loc_bootstrap_marker_1 = numeric(B)
        loc_bootstrap_marker_2 = numeric(B)
        kernel_opt_bootstrap_marker_1 = numeric(B)
        kernel_opt_bootstrap_marker_2 = numeric(B)
        kernel_hscv_bootstrap_marker_1 = numeric(B)
        kernel_hscv_bootstrap_marker_2 = numeric(B)
        
        for (b in 1:B) {
          indexes_sanos = sample(1:n, n, replace = TRUE)
          indexes_enfermos = sample(1:n, n, replace = TRUE)
          
          sanos_marker_1_b = sanos_marker_1[indexes_sanos]
          sanos_marker_2_b = sanos_marker_2[indexes_sanos]
          
          enfermos_marker_1_b = enfermos_marker_1[indexes_enfermos]
          enfermos_marker_2_b = enfermos_marker_2[indexes_enfermos]
          
          # AUC
          auc_bootstrap_marker_1[b] = parametric_auc(sanos_marker_1_b, enfermos_marker_1_b, bc = boxcox)
          auc_bootstrap_marker_2[b] = parametric_auc(sanos_marker_2_b, enfermos_marker_2_b, bc = boxcox)
          
          # Youden
          youden_bootstrap_marker_1[b] = parametric_youden(sanos_marker_1_b, enfermos_marker_1_b, bc = boxcox)
          youden_bootstrap_marker_2[b] = parametric_youden(sanos_marker_2_b, enfermos_marker_2_b, bc = boxcox)
          
          # Parametric
          binormal_bootstrap_marker_1[b] = EtaBinormal(sanos_marker_1_b, enfermos_marker_1_b, bc = boxcox)
          binormal_bootstrap_marker_2[b] = EtaBinormal(sanos_marker_2_b, enfermos_marker_2_b, bc = boxcox)
          
          # LoC
          loc_bootstrap_marker_1[b] = parametric_loc(sanos_marker_1_b, enfermos_marker_1_b, bc = boxcox)
          loc_bootstrap_marker_2[b] = parametric_loc(sanos_marker_2_b, enfermos_marker_2_b, bc = boxcox)
          
          # Kernel opt
          kernel_opt_bootstrap_marker_1[b] = EtaKernel(sanos_marker_1_b, enfermos_marker_1_b, "optimo", mesh_size = 300, bc = boxcox)
          kernel_opt_bootstrap_marker_2[b] = EtaKernel(sanos_marker_2_b, enfermos_marker_2_b, "optimo", mesh_size = 300, bc = boxcox)
          
          # Kernel hscv
          kernel_hscv_bootstrap_marker_1[b] = EtaKernel(sanos_marker_1_b, enfermos_marker_1_b, "hscv", mesh_size = 300, bc = boxcox)
          kernel_hscv_bootstrap_marker_2[b] = EtaKernel(sanos_marker_2_b, enfermos_marker_2_b, "hscv", mesh_size = 300, bc = boxcox)
        }
        
        # Cálculo de diferencias entre los dos marcadores
        delta_auc = auc_bootstrap_marker_2 - auc_bootstrap_marker_1
        delta_youden = youden_bootstrap_marker_2 - youden_bootstrap_marker_1
        delta_parametric = binormal_bootstrap_marker_2 - binormal_bootstrap_marker_1
        delta_loc = loc_bootstrap_marker_2 - loc_bootstrap_marker_1
        delta_kernel_opt = kernel_opt_bootstrap_marker_2 - kernel_opt_bootstrap_marker_1
        delta_kernel_hscv = kernel_hscv_bootstrap_marker_2 - kernel_hscv_bootstrap_marker_1
        
        p_val_auc[i] = 2 * min(mean(delta_auc < 0), mean(delta_auc > 0))
        p_val_youden[i] = 2 * min(mean(delta_youden < 0), mean(delta_youden > 0))
        p_val_param[i] = 2 * min(mean(delta_parametric < 0), mean(delta_parametric > 0))
        p_val_loc[i] = 2 * min(mean(delta_loc < 0), mean(delta_loc > 0))
        p_val_kernel_opt[i] = 2 * min(mean(delta_kernel_opt < 0), mean(delta_kernel_opt > 0))
        p_val_kernel_hscv[i] = 2 * min(mean(delta_kernel_hscv < 0), mean(delta_kernel_hscv > 0))
      }
      
      potencia_auc = mean(p_val_auc < alpha)
      potencia_youden = mean(p_val_youden < alpha)
      potencia_param = mean(p_val_param < alpha)
      potencia_loc = mean(p_val_loc < alpha)
      potencia_kernel_opt = mean(p_val_kernel_opt < alpha)
      potencia_kernel_hscv = mean(p_val_kernel_hscv < alpha)
      
      covariance_results[[paste0("N_", n)]] = list(
        auc = potencia_auc,
        youden = potencia_youden,
        param = potencia_param,
        loc_p = potencia_loc,
        kernel_opt = potencia_kernel_opt,
        kernel_hscv = potencia_kernel_hscv
      )
    }
    
    results[[paste0("correlation_", covariance_enfermos)]] = covariance_results
  }
  
  json = toJSON(results, pretty = T, digits = NA)
  dir_path = here('two_biomarkers', 'jsons')
  filename = paste0('potencias_', num,'.json')
  full_path = file.path(dir_path, filename )
  write(json, file = full_path)

}

scenario_n_two_biomarkers(0.36, 0.56, 1, 4, FALSE, 1)
scenario_n_two_biomarkers(0.74, 1.17, 1, 4, FALSE, 2)
scenario_n_two_biomarkers(1.19, 1.88, 1, 4, FALSE, 3)
scenario_n_two_biomarkers(1.81, 2.87, 1, 4, FALSE, 4)


scenario_n_two_biomarkers(0.36, 0.74, 1, 1, FALSE, 5)
scenario_n_two_biomarkers(0.36, 1.19, 1, 1, FALSE, 6)
scenario_n_two_biomarkers(0.36, 1.81, 1, 1, FALSE, 7)
scenario_n_two_biomarkers(0.74, 1.19, 1, 1, FALSE, 8)
scenario_n_two_biomarkers(0.74, 1.81, 1, 1, FALSE, 8)
scenario_n_two_biomarkers(1.19, 1.81, 1, 1, FALSE, 8)


