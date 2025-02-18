import os
import json
os.chdir("./niveles_y_potencias/Scenarios/potencias/jsons")

file_groups = {
    "tabla1": ["potencias_1.json", "potencias_2.json", "potencias_3.json"],
    "tabla2": ["potencias_4.json", "potencias_5.json", "potencias_6.json", "potencias_7.json"],
    "tabla3": ["potencias_8.json", "potencias_9.json", "potencias_10.json"]
}


placeholder = "."


columns = ["20", "50", "100"]
row_labels = [
    "$\\widehat{AUC}$",
    "$\\widehat{J}$",
    "$\\widehat{LoC}$",
    "$\\hat{\\eta}_{log}^{N}$",
    "$\\hat{\\eta}_{log}^{K = g, h=csv}$",
    "$\\hat{\\eta}_{log}^{K = g, h = h*}$"
]


metric_to_json_key = {
    "$\\widehat{AUC}$": "auc",
    "$\\widehat{J}$": "youden",
    "$\\hat{\\eta}_{log}^{NonParametric}$": "no_param",
    "$\\widehat{LoC}$": "LoC",
    "$\\hat{\\eta}_{log}^{N}$": "param",
    "$\\hat{\\eta}_{log}^{K = g, h=csv}$": "kernel_hscv",
    "$\\hat{\\eta}_{log}^{K = g, h = h*}$": "kernel_opt"
}


def read_json(file_path):
    if not os.path.exists(file_path):
        return {f"N_{n}": {k: placeholder for k in metric_to_json_key.values()} for n in [20, 50, 100]}
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def generate_latex_table(group_files, distributions):
    latex_content = "\\begin{table}[H]\n\\centering\n"
    latex_content += "\\caption{Resultados de potencias para diferentes estimadores.}\n"
    latex_content += "\\renewcommand{\\arraystretch}{1.5} % Incrementa el espaciado entre filas\n"
    latex_content += "\\begin{tabular}{lllccc}\n"
    latex_content += "\\multicolumn{3}{c}{}  & \\multicolumn{3}{c}{\\(n_x = n_y\\)} \\\\\n"
    latex_content += "\\cmidrule(lr){4-6}\n"
    latex_content += "Distribución sanos & Distribución enfermos &   & 20 & 50 & 100 \\\\\n"
    latex_content += "\\midrule\n"
    
    for dist_sanos, dist_enfermos, file_name in distributions:
        json_data = read_json(file_name)
        for row_label in row_labels:
            row = f" {dist_sanos if row_label == '$\\widehat{AUC}$' else ''} & "
            row += f"{dist_enfermos if row_label == '$\\widehat{AUC}$' else ''} & "
            row += f"{row_label} & "

            row += " & ".join([str(json_data.get(f"N_{n}", {}).get(metric_to_json_key[row_label], placeholder)) for n in [20, 50, 100]])
            row += " \\\\\n"
            latex_content += row
        latex_content += "&&&&&\\\\\n"  # Espaciado entre distribuciones

    latex_content += "\\bottomrule\n"
    latex_content += "\\end{tabular}\n"
    latex_content += "\\label{tab:simulated_power}\n"
    latex_content += "\\end{table}\n"
    return latex_content


distributions = {
    "tabla1": [
        ("N(0,1)", "N(0.36,1)", file_groups["tabla1"][0]),
        ("", "N(0.44,1.4)", file_groups["tabla1"][1]),
        ("", "N(0.26,0.3)", file_groups["tabla1"][2])
    ],
    "tabla2": [
        ("LN(0,1)", "LN(0.28,0.5)", file_groups["tabla2"][0]),
        ("", "LN(0.46,3/2)", file_groups["tabla2"][1]),
        ("", "LN(0.26,0.2)", file_groups["tabla2"][2]),
        ("", "LN(0.57,2)", file_groups["tabla2"][3])
    ],
    "tabla3": [
        ("Gamma(0.5,0.5)", "Gamma(1.05,1)", file_groups["tabla3"][0]),
        ("", "Gamma(3.28,4)", file_groups["tabla3"][1]),
        ("", "Gamma(0.35,1/8)", file_groups["tabla3"][2])
    ]
}

# Generar y guardar las tablas en un archivo LaTeX
output_file = "tablas_latex_loc.txt"
with open(output_file, "w", encoding="utf-8") as f:
    for table_name, table_distributions in distributions.items():
        latex_table = generate_latex_table(file_groups[table_name], table_distributions)
        f.write(latex_table + "\n\n")

print(f"Las tablas se han guardado en el archivo {output_file}.")
