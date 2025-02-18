import json
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

os.chdir("./niveles_y_potencias/two_biomarkers/potencias/jsons")

def json_to_csv(file_groups, output_csv):
    metric_names = {
        "youden": "Youden",
        "loc_p": "LoC",
        "param": "Eta (paremétrico)",
        "kernel_hscv": "Eta (kernel hscv)",
        "auc": "AUC",
        "kernel_opt": "Eta (kernel h*)"
    }
    
    rows = []
    metrics_interes = list(metric_names.keys())
    nombres_grupos = ["cruzan", "no cruzan"]
    
    for group_idx, file_list in enumerate(file_groups):
        group_name = f"Curvas ROC: {nombres_grupos[group_idx]}"
        for file in file_list:
            with open(file, 'r') as f:
                data = json.load(f)

                for key, corr_data in data.items():
                    if key.startswith("correlation_"):

                        correlation_val = key.split("_")[1]

                        for n_key, metrics in corr_data.items():
                            if n_key.startswith("N_"):

                                N_val = n_key.split("_")[1]

                                for metric, value in metrics.items():
                                    if metric in metrics_interes:
                                        rows.append({
                                            "group": group_name,
                                            "correlation": correlation_val,
                                            "N": N_val,
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
    plt.title("Potencias del contraste con T", fontsize=14)
    plt.ylabel("Potencia", fontsize=12)
    plt.xlabel("")
    plt.xticks(rotation=45)
    

    ax.legend(title="Grupo", fontsize=10, title_fontsize=12,
              bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.subplots_adjust(right=0.75)
    plt.show()


crossing    = ["potencias_1.json", "potencias_2.json", "potencias_3.json", "potencias_4.json"] 
superior    = ["potencias_5.json", "potencias_6.json", "potencias_7.json", "potencias_8.json", "potencias_9.json", "potencias_10.json"] 

df = json_to_csv([crossing, superior], "resultados.csv")
plot_boxplot("resultados.csv")
