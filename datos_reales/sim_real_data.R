library(jsonlite)
library(readr) 
library(dplyr)
library(here)

sim_bootstrap = function(m1_sanos, m1_enfermos, m2_sanos, m2_enfermos, name) {
  nsanos = length(m1_sanos)
  nenfermos = length(m1_enfermos)
  B = 10000
  
  auc_1 = numeric(B)
  auc_2 = numeric(B)
  
  j_1 = numeric(B)
  j_2 = numeric(B)
  
  loc_1 = numeric(B)
  loc_2 = numeric(B)
  
  eta_param_1 = numeric(B)
  eta_param_2 = numeric(B)
  
  eta_kernel_opt_1 = numeric(B)
  eta_kernel_opt_2 = numeric(B)
  
  eta_kernel_hscv_1 = numeric(B)
  eta_kernel_hscv_2 = numeric(B)
  
  for (b in 1:B) {
    indexes_sanos = sample(1:nsanos, nsanos, replace = TRUE)
    indexes_enfermos = sample(1:nenfermos, nenfermos, replace = TRUE)
    
    sanos_1_b = m1_sanos[indexes_sanos]
    sanos_2_b = m2_sanos[indexes_sanos]
    
    enfermos_1_b = m1_enfermos[indexes_enfermos]
    enfermos_2_b = m2_enfermos[indexes_enfermos]
    
    auc_1[b] = parametric_auc(sanos_1_b, enfermos_1_b, bc = TRUE)
    auc_2[b] = parametric_auc(sanos_2_b, enfermos_2_b, bc = TRUE)
    
    j_1[b] = parametric_youden(sanos_1_b, enfermos_1_b, bc = TRUE)
    j_2[b] = parametric_youden(sanos_2_b, enfermos_2_b, bc = TRUE)
    
    loc_1[b] = parametric_loc(sanos_1_b, enfermos_1_b, bc = TRUE)
    loc_2[b] = parametric_loc(sanos_2_b, enfermos_2_b, bc = TRUE)
    
    eta_param_1[b] = EtaBinormal(sanos_1_b, enfermos_1_b, bc = TRUE)
    eta_param_2[b] = EtaBinormal(sanos_2_b, enfermos_2_b, bc = TRUE)
    
    eta_kernel_opt_1[b] = EtaKernel(sanos_1_b, enfermos_1_b, "optimo", mesh_size = 300, bc = TRUE)
    eta_kernel_opt_2[b] = EtaKernel(sanos_2_b, enfermos_2_b, "optimo", mesh_size = 300, bc = TRUE)
    
    eta_kernel_hscv_1[b] = EtaKernel(sanos_1_b, enfermos_1_b, "hscv", mesh_size = 300, bc = TRUE)
    eta_kernel_hscv_2[b] = EtaKernel(sanos_2_b, enfermos_2_b, "hscv", mesh_size = 300, bc = TRUE)
  }
  delta_auc = auc_2 - auc_1
  delta_j = j_2 - j_1
  delta_loc = loc_2 - loc_1
  delta_param = eta_param_2 - eta_param_1
  delta_k_opt = eta_kernel_opt_2 - eta_kernel_opt_1
  delta_k_hscv = eta_kernel_hscv_2 - eta_kernel_hscv_1
  
  p_val_auc = 2 * min(mean(delta_auc < 0), mean(delta_auc > 0))
  p_val_youden = 2 * min(mean(delta_j < 0), mean(delta_j > 0))
  p_val_loc = 2 * min(mean(delta_loc < 0), mean(delta_loc > 0))
  p_val_param = 2 * min(mean(delta_param < 0), mean(delta_param > 0))
  p_val_kernel_opt = 2 * min(mean(delta_k_opt < 0), mean(delta_k_opt > 0))
  p_val_kernel_hscv = 2 * min(mean(delta_k_hscv < 0), mean(delta_k_hscv > 0))
  
  res = list(
    p_auc = p_val_auc,
    p_youden = p_val_youden,
    p_param = p_val_param,
    p_loc = p_val_loc,
    p_opt = p_val_kernel_opt,
    p_hscv = p_val_kernel_hscv
  )
  
  output_dir = here("datos_reales", "jsons")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  output_path = file.path(output_dir, paste0(name, ".json"))
  write_json(res, output_path, pretty = TRUE)
}

# Importar los datos del primer archivo
data = read_csv(here("datos_reales", "final_data_1.csv"))
ca125_sanos = data$ca125[data$status == 0]
ca125_enfermos = data$ca125[data$status == 1]
ca199_sanos = data$ca199[data$status == 0]
ca199_enfermos = data$ca199[data$status == 1]

# Importar los datos del segundo archivo
data = read_csv(here("datos_reales", "final_data_2.csv"))
x_sanos = data$x_value[data$x_status == 0]
x_enfermos = data$x_value[data$x_status == 1]
y_sanos = data$y_value[data$y_status == 0]
y_enfermos = data$y_value[data$y_status == 1]

# Importar los datos del tercer archivo
data = read_csv(here("datos_reales", "zhang.csv"))
x_sanos_z = as.numeric(data$modality1[data$status...4 == 0])
x_enfermos_z = as.numeric(data$modality1[data$status...4 == 1])
y_sanos_z = as.numeric(data$modality2[data$status...4 == 0])
y_enfermos_z = as.numeric(data$modality2[data$status...4 == 1])

# Ejecutar la simulación con los datos cargados
sim_bootstrap(ca125_sanos, ca125_enfermos, ca199_sanos, ca199_enfermos, "data_1")
sim_bootstrap(x_sanos, x_enfermos, y_sanos, y_enfermos, "data_2")
sim_bootstrap(x_sanos_z, x_enfermos_z, y_sanos_z, y_enfermos_z, "data_zhang")

