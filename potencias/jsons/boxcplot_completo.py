import json
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir("./niveles_y_potencias/Scenarios/potencias/jsons")

def json_to_csv(file_groups, output_csv):
    """
    Procesa una lista de listas de archivos JSON.
    Cada sublista se tratará como un grupo distinto y se añadirá una columna "group" al dataframe.
    """
    metric_names = {
        "youden": "Youden",
        "LoC": "LoC",
        "param": "Eta (paremétrico)",
        "kernel_hscv": "Eta (kernel hscv)",
        "auc": "AUC",
        "kernel_opt": "Eta (kernel h*)"
    }
    rows = []
    for group_idx, file_list in enumerate(file_groups):
        nombres = ["normales", "lognormales", "gamma"]
        group_name = f"Poblaciones {nombres[group_idx]}"
        for file in file_list:
            with open(file, 'r') as f:
                data = json.load(f)
                for key, values in data.items():
                    if key.startswith("N_"):
                        for metric, value in values.items():
                            if metric in ["youden", "param", "kernel_hscv", "auc", "kernel_opt", "LoC"]:
                                rows.append({
                                    "group": group_name,
                                    "metric": metric_names.get(metric, metric),
                                    "value": value
                                })
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    return df

def plot_boxplot(csv_file):
    df = pd.read_csv(csv_file)
    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid")
    
    ax = sns.boxplot(x="metric", y="value", hue="group", data=df, palette="Blues")
    plt.title("Potencias agrupadas por distribución y medida resumen", fontsize=14)
    plt.ylabel("Potencia", fontsize=12)
    plt.xlabel("")
    plt.xticks(rotation=45)
    
    ax.legend(title="Grupo", fontsize=10, title_fontsize=12,
              bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    
    plt.subplots_adjust(right=0.75)
    plt.show()



normales    = ["potencias_1.json", "potencias_2.json", "potencias_3.json"] 
lognormales = ["potencias_4.json", "potencias_5.json", "potencias_6.json", "potencias_7.json"]
gammas      = ["potencias_8.json", "potencias_9.json", "potencias_10.json"] 

df = json_to_csv([normales, lognormales, gammas], "resultados.csv")
plot_boxplot("resultados.csv")

