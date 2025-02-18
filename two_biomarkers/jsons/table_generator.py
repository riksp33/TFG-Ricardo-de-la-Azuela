import json
import os

# Cambia el directorio al lugar donde se encuentran los JSON
os.chdir('./niveles_y_potencias/two_biomarkers/potencias')

def json_to_latex_longtable(json_files, output_file, auc_values, metrics):
    # AUC y Métrica van en pares
    measures = [
        "$\widehat{AUC}$",
        "$\widehat{J}$",
        "$\widehat{LoC}$",
        "$\hat{\eta}_{log}^{N}$",
        "$\hat{\eta}_{log}^{K = g, h=csv}$",
        "$\hat{\eta}_{log}^{K = g, h = h*}$"
    ]

    # Inicializa el código LaTeX
    latex_code = []
    latex_code.append("\\renewcommand{\\arraystretch}{1.5}")
    latex_code.append("\\begin{longtable}{lllcccccc}")
    latex_code.append("\caption{Potencias para marcadores bottom} \label{tab:simulated_power3} \\\\")
    latex_code.append("\\toprule")
    latex_code.append("\multicolumn{3}{c}{}  & \multicolumn{6}{c}{\(n_x = n_y\)} \\\\")
    latex_code.append("\cmidrule(lr){4-9}")
    latex_code.append(" &  &  & \multicolumn{3}{c}{$\\rho = 0.0$} & \multicolumn{3}{c}{$\\rho = 0.5$} \\\\")
    latex_code.append("\cmidrule(lr){4-6} \cmidrule(lr){7-9}")
    latex_code.append(" A  &    B      &        & 20 & 50 & 100 & 20 & 50 & 100 \\\\")
    latex_code.append("\midrule")
    latex_code.append("\endfirsthead")

    latex_code.append("\caption[]{Continúa desde la página anterior} \\\\")
    latex_code.append("\\toprule")
    latex_code.append("\multicolumn{3}{c}{}  & \multicolumn{6}{c}{\(n_x = n_y\)} \\\\")
    latex_code.append("\cmidrule(lr){4-9}")
    latex_code.append(" &  &  & \multicolumn{3}{c}{$\\rho = 0.0$} & \multicolumn{3}{c}{$\\rho = 0.5$} \\\\")
    latex_code.append("\cmidrule(lr){4-6} \cmidrule(lr){7-9}")
    latex_code.append(" A  &    B      &        & 20 & 50 & 100 & 20 & 50 & 100 \\\\")
    latex_code.append("\midrule")
    latex_code.append("\endhead")

    latex_code.append("\\bottomrule")
    latex_code.append("\endfoot")

    # Itera sobre los archivos JSON y genera filas
    for auc_idx, json_file in enumerate(json_files):
        with open(json_file, 'r') as file:
            data = json.load(file)

        auc = auc_values[auc_idx]
        metric = metrics[auc_idx]

        first_row = True  # Para manejar la primera fila del bloque
        for measure_idx, measure in enumerate(measures):
            row_values_corr_0 = []
            row_values_corr_5 = []
            for sample_size in ["20", "50", "100"]:
                sample_key = f"N_{sample_size}"

                # Procesa correlación 0.0
                if "correlation_0.0" in data:
                    corr_0_data = data["correlation_0.0"].get(sample_key, {})
                    if measure_idx == 0:
                        row_values_corr_0.append(corr_0_data.get("auc", "N/A"))
                    elif measure_idx == 1:
                        row_values_corr_0.append(corr_0_data.get("youden", "N/A"))
                    elif measure_idx == 2:
                        row_values_corr_0.append(corr_0_data.get("loc_p", "N/A"))
                    elif measure_idx == 3:
                        row_values_corr_0.append(corr_0_data.get("param", "N/A"))
                    elif measure_idx == 4:
                        row_values_corr_0.append(corr_0_data.get("kernel_hscv", "N/A"))
                    elif measure_idx == 5:
                        row_values_corr_0.append(corr_0_data.get("kernel_opt", "N/A"))

                # Procesa correlación 0.5
                if "correlation_0.5" in data:
                    corr_5_data = data["correlation_0.5"].get(sample_key, {})
                    if measure_idx == 0:
                        row_values_corr_5.append(corr_5_data.get("auc", "N/A"))
                    elif measure_idx == 1:
                        row_values_corr_5.append(corr_5_data.get("youden", "N/A"))
                    elif measure_idx == 2:
                        row_values_corr_5.append(corr_5_data.get("loc_p", "N/A"))
                    elif measure_idx == 3:
                        row_values_corr_5.append(corr_5_data.get("param", "N/A"))
                    elif measure_idx == 4:
                        row_values_corr_5.append(corr_5_data.get("kernel_hscv", "N/A"))
                    elif measure_idx == 5:
                        row_values_corr_5.append(corr_5_data.get("kernel_opt", "N/A"))

            # Formato LaTeX para la fila
            if first_row:
                latex_code.append(
                    f"{auc} & {metric} & {measure} & "
                    f"{' & '.join(map(str, row_values_corr_0))} & {' & '.join(map(str, row_values_corr_5))} \\\\"
                )
                first_row = False
            else:
                latex_code.append(
                    f"    &       & {measure} & "
                    f"{' & '.join(map(str, row_values_corr_0))} & {' & '.join(map(str, row_values_corr_5))} \\\\"
                )
        latex_code.append("\\\\[0.2cm]")  # Espaciado vertical entre bloques

    latex_code.append("\end{longtable}")

    # Guarda el resultado en el archivo de salida
    with open(output_file, 'w') as file:
        file.write('\n'.join(latex_code))

# Ejemplo de uso
tabla_0 = ["potencias_1.json", "potencias_2.json", "potencias_3.json", "potencias_4.json"]
tabla_1 = ["potencias_5.json", "potencias_6.json", "potencias_7.json", "potencias_8.json", "potencias_9.json", "potencias_10.json"]
auc_values = ["0.6", "0.6", "0.6", "0.7", "0.7", "0.8"]
metrics = ["0.7", "0.8", "0.9", "0.8", "0.9", "0.9"]
json_files = ['./jsons/' + file for file in tabla_1]  # Ajusta la ruta según sea necesario
output_file = "tabla_latex_bottom.txt"
json_to_latex_longtable(json_files, output_file, auc_values, metrics)

json_files = ['./jsons/' + file for file in tabla_0]
output_file = "tabla_latex_top.txt"
json_to_latex_longtable(json_files, output_file, auc_values, metrics)
