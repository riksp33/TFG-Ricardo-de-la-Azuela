import json
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir("./bias_rmse/jsons")

def json_to_csv(files, output_csv):
    """
    Procesa una lista de archivos JSON con el siguiente formato:
    
    {
      "header": [...],
      "AUC:0.6": {
         "size:20": {
            "eta_pob": [...],
            "no_param": { ... },
            "param": {
              "bias": [valor],
              "rmse": [valor]
            },
            "kernel_opt": {
              "bias": [valor],
              "rmse": [valor]
            },
            "kernel_hscv": {
              "bias": [valor],
              "rmse": [valor]
            }
         },
         "size:50": { ... },
         "size:100": { ... }
      },
      "AUC:0.75": { ... },
      "AUC:0.9": { ... }
    }
    
    Para cada archivo, se extraen los valores de _bias_ y _rmse_ para los métodos
    "kernel_hscv", "kernel_opt" y "param". Se guarda un registro por cada valor.
    """
    methods_of_interest = ["kernel_hscv", "kernel_opt", "param"]
    method_names = {
        "kernel_hscv": "Eta (kernel hscv)",
        "kernel_opt": "Eta (Kernel h*)",
        "param": "Eta (paramétrico)"
    }
    measure_names = {
        "bias": "Bias",
        "rmse": "RMSE"
    }
    
    rows = []
    for file in files:
        with open(file, 'r') as f:
            data = json.load(f)
            for auc_key, auc_data in data.items():
                if auc_key.startswith("AUC:"):
                    auc_val = auc_key.split(":")[1] 
                    for size_key, size_data in auc_data.items():
                        if size_key.startswith("size:"):
                            size_val = size_key.split(":")[1]
                            for method in methods_of_interest:
                                if method in size_data:
                                    method_data = size_data[method]
                                    for measure in ["bias", "rmse"]:
                                        if measure in method_data:
                                            for value in method_data[measure]:
                                                rows.append({
                                                    "AUC": auc_val,
                                                    "size": size_val,
                                                    "method": method_names.get(method, method),
                                                    "measure": measure_names.get(measure, measure),
                                                    "value": value
                                                })
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    return df

def plot_boxplot(csv_file):
    """
    Lee el CSV generado y crea un gráfico de box plots en el que se representan,
    para cada método (Kernel HSCV, Kernel óptimo, Parámetro), dos cajas correspondientes a Bias y RMSE.
    """
    df = pd.read_csv(csv_file)
    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")
    
    ax = sns.boxplot(x="method", y="value", hue="measure", data=df, palette="Blues")
    
    plt.title("Bias y RMSE para diferentes métodos. Poblaciones normales", fontsize=14)
    plt.ylabel("", fontsize=12)
    plt.xlabel("", fontsize=12)
    plt.xticks(rotation=45)
    
    ax.legend(title="Medida", fontsize=10, title_fontsize=12,
              bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.subplots_adjust(right=0.75)
    plt.show()

files = [f"tabla{i}" for i in range(1, 11)]

df = json_to_csv(files, "resultados.csv")
plot_boxplot("resultados.csv")
