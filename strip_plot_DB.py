# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
import sys
import os
import re
import argparse

# FUNCION PARA SELECCIONAR LOS ARCHIVOS DE LOS ESTUDIOS
def select_archives():
    files = []
    # El segundo argumento del script es la ruta de la carpeta que contiene las matrices de expresión diferencial
    # y el tercer argumento es la ruta del directorio de salida
    parser = argparse.ArgumentParser(description = "Procesa un archivo de entrada")
    parser.add_argument("-i", "--input", required = True, help = "Ruta del archivo de entrada")
    parser.add_argument("-o", "--output", required = True, help = "Ruta del archivo de salida")
    args = parser.parse_args()
    file_path = args.input
    output_path = args.output
    # Normaliza el path para evitar problemas con rutas relativas o absolutas
    file_path_clean = os.path.normpath(os.path.abspath(file_path))
    output_path = os.path.normpath(os.path.abspath(output_path))
    # Busca todas las carpetas cuyo nombre empiece por "SRP" y, dentro de cada una,
    # añade todas las .tsv que empiecen por DESeq2, edgeR o limma.
    metodos_prefijos = ("DESeq2", "edgeR", "limma")
    contador_carpetas = 1
    for root, dirs, files_in_dir in os.walk(file_path_clean):
        # Identifica directorios de estudio (ej: SRP123456)
        if os.path.basename(root).startswith("SRP"):
            # Recorre las subcarpetas del estudio en busca de la carpeta de regresion
            for subroot, subdirs, subfiles in os.walk(root):
                    # Añade solo los archivos .tsv cuyos nombres empiecen por alguno de los metodos y en los que se enfrenten 2 condiciones (vs)
                    for file in subfiles:
                        if file.endswith(".tsv") and file.startswith(metodos_prefijos) and ("vs" in file):
                            file_path_full = os.path.join(subroot, file)
                            files.append(file_path_full)
    return files, output_path

# FUNCIÓN PARA CREAR EL DATAFRAME CON LA INFORMACIÓN DE LOS ARCHIVOS SELECCIONADOS
def data_frame(files, percentil):
    # Lista con los métodos
    metodos = ("DESeq2", "edgeR", "limma")
    # Lista con las bases de datos
    dbs = ["miRBase", "mirgenedb", "mircarta"]
    # Tamaño de cada base de datos (número de miRNAs que contienen) para la especie humana
    db_size = count_mirnas_DB()

    # Lee cada archivo seleccionado, cuenta el número de miRNAs significativos (padj < 0.05)... Y guarda esta información en un DataFrame
    # que contiene el nombre del estudio (folder_name), la base de datos usada, el método (ambos extraidos de files_names),
    #  el número de miRNAs significativos encontrados en cada archivo y otra información relevante 
    data = []
    for file in files:
        # Primero se lee el archivo en un DataFrame y se obtienen valores iniciales de interés
        df = pd.read_csv(file, sep = "\t")
        total_mirnas = len(df)
        significant_mirnas = df[df["padj"] < 0.05]
        significant_mirnas_log2fc = df[(df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)]
        num_significant_mirnas = len(significant_mirnas)
        num_significant_mirnas_log2fc = len(significant_mirnas_log2fc)

        # Extrae la base de datos y el método del nombre del archivo (asumiendo que el formato es "metodo_TEXTO_DB.txt")
        metodo = "UNDEFINED"
        base_de_datos = "UNDEFINED"
        # Busca el método, la base de datos y el nombre del estudio en la ruta del archivo
        ruta = file.split(os.sep)
        for part in ruta:
            if part.startswith("SRP"):
                # Si empieza con SRP se asume que es el nombre del estudio, sigue la estructura "SRP_DB"
                estudio_db = part.split("_")
                estudio = estudio_db[0]
                if estudio_db[1] in dbs:
                    base_de_datos = estudio_db[1]
            if part.startswith(metodos):
                # Si empieza con método se asume que es el método, sigue la estructura "METODO_....tsv"
                metodo_etc = part.split("_")
                metodo = metodo_etc[0]
            if part.startswith("de_"):
                # Si empieza con de_ se asume que es el metodo de asignación de lecturas, sigue la estructura "de_rcadj" o "de_rcsa"
                assignment = part

        # Para el cálculo de porcentajes se utilizará num_significant_mirnas_log2fc ya que son los que tienen una expresión diferencial y significativa
        if total_mirnas > 0:
            estudio_percent = num_significant_mirnas_log2fc / total_mirnas * 100
            db_percent = num_significant_mirnas_log2fc / db_size[base_de_datos] * 100
        else:
            estudio_percent = None
            db_percent = None

        # Calcula el percentil del valor absoluto de log2FoldChange para los miRNAs significativos
        percentil_abs_log2fc = significant_mirnas["log2FoldChange"].abs().quantile(percentil / 100.0)

        # Agrega la información a data como una lista de diccionarios
        data.append({"estudio": estudio, "db": base_de_datos, "metodo": metodo, "assignment": assignment, "miRNAs_significativos": num_significant_mirnas, "miRNAs_significativos_DE": num_significant_mirnas_log2fc, "total_miRNAs_testeados": total_mirnas, "porcentaje_sig_estudio": estudio_percent, "porcentaje_sig_db": db_percent, "abs(FoldChange)_Percentil_sig": percentil_abs_log2fc})
    # Crea un DataFrame a partir de la lista de diccionarios
    df = pd.DataFrame(data)
    df_sorted = df.sort_values(by = ["estudio", "db", "metodo"]).reset_index(drop = True)
    # Guarda el DataFrame en un archivo TSV para verificar que se ha creado correctamente
    df_sorted.to_csv(f"data_frame_strip_plot.tsv", sep = "\t", index = False)  
    return df_sorted

# FUNCION PARA CONTAR EL NÚMERO DE MIRNAS EN LAS BASES DE DATOS
def count_mirnas_DB(files = "/shared/bak/TFG/serrano/db", especie = "Hsa"):
    # Cuenta el número de miRNAs en cada base de datos para la especie dada (hsa)
    # Hay que tener en cuenta que cuenta el número de miRNAs maduros y no el preName
    db_counts = {}
    count = 0
    for root, dirs, files_in_dir in os.walk(files):
        for file in files_in_dir:
            if file.endswith(".tsv"):
                file_path_full = os.path.join(root, file)
                df = pd.read_csv(file_path_full, sep="\t")
                # Se asume que el formato del archivo es "base_de_datos.tsv" y que la base de datos es el nombre del archivo sin la extensión
                base_de_datos = os.path.splitext(file)[0]
                # En la columna "matureNameString" se cuenta el número de ocurrencias de "especie" (puede aparecer varias veces en la misma celda)
                # Usamos una búsqueda insensible a mayúsculas/minúsculas y escapamos la cadena para evitar metacaracteres
                especie = re.escape(especie)
                # Se usa "fr"(?i){pattern}"" para hacer la búsqueda insensible a mayúsculas/minúsculas
                count = df["matureNameString"].str.count(fr"(?i){especie}").sum()
                db_counts[base_de_datos] = count
    return db_counts

# FUNCIÓN UNIFICADA PARA CREAR STRIP PLOTS LAS DISTINTAS BASES DE DATOS, MÉTODOS Y ESTUDIOS, CON LOS MIARNAS SIGNIFICATIVOS
def strip_plot_padj(df, type):
    # Primero se crea una copia del DataFrame original y se representa el strip plot con seaborn
    df_plot = df.copy()

    # Generamos los subplots
    fig, axex = plt.subplots(1, 2, figsize = (10, 8), sharey = True)

    # Se establecen los colores para los boxplots y los puntos de cada método, y se configura el estilo de la gráfica con seaborn
    sns.set_theme(style = "ticks")
    box_colors = {"DESeq2": "#9ecae1", "limma": "#a1d99b", "edgeR": "#fdae6b"}
    point_colors = {"DESeq2": "#3182bd", "limma": "#31a354", "edgeR": "#e6550d"}

    # Se obtiene la lista de asignaciones únicas para iterar sobre ellas
    asignacion = df_plot["assignment"].unique().tolist()

    # Dependiendo de el "type" que se elija se representa porcentaje_sig_estudio o porcentaje_sig_db
    for i, asing in enumerate(asignacion):
        ax = axex[i]

        # Pintamos el boxplot de fondo para cada método, con un color específico para cada uno
        sns.boxplot(data = df_plot[df_plot["assignment"] == asing], x = "db", y = type, hue = "metodo", palette = box_colors, width = 0.6, dodge = True, fliersize = 0, boxprops = dict(alpha = 0.6), ax = ax)

        # Y pintamos los puntos por encima, con un color específico para cada método
        sns.stripplot(data = df_plot[df_plot["assignment"] == asing], x = "db", y = type, hue = "metodo", palette = point_colors, dodge = True, jitter = 0.15, size = 5, alpha = 0.8, marker = "o", ax = ax)

        # Se añade el titulo a cada subplot y se eliminan marcas individuales de los ejes 
        ax.set_title(f"{asing}", fontweight = "bold", pad = 12)
        ax.legend_.remove()
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Nombramos cada eje y la gráfica, añadiendo también una leyenda para los métodos
    fig.supxlabel("Base de datos")
    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles[:3], labels[:3], loc = "upper left", ncol = 1, title = "Marco de Expresion Diferencial", title_fontsize = "11", fontsize = "9", bbox_to_anchor = (1.02, 1), frameon = True)

    if "db" in type:
        fig.suptitle("Porcentaje miRNAs significativos diferencialmente expresados respecto a la base de datos", fontweight = "bold")
        fig.supylabel("Porcentaje de miRNAs respecto a la base de datos")
        sns.despine()
        plt.tight_layout()
        plt.savefig("porcentaje_respecto_db.png", dpi = 300, bbox_inches = "tight")
    elif "DE" in type:
        fig.suptitle("Total miRNAs significativos diferencialmente expresados", fontweight = "bold")
        fig.supylabel("Número de miRNAs")
        sns.despine()
        plt.tight_layout()
        plt.savefig("total_mirnas_DE.png", dpi = 300, bbox_inches = "tight")
    else:
        fig.suptitle("Porcentaje miRNAs significativos diferencialmente expresados respecto al total evaluados", fontweight = "bold")
        fig.supylabel("Porcentaje de miRNAs respecto total de miRNAs analizado")
        sns.despine()
        plt.tight_layout()
        plt.savefig("porcentaje_respecto_testeados.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

# FUNCIÓN PARA CREAR STRIP PLOTS DE LAS DISTINTAS BASES DE DATOS, MÉTODOS Y ESTUDIOS, USANDO EL PERCENTIL 80 DE LOS VALORES ABSOLUTOS DE LOG2FC
def strip_plot_percentil_log2fc(df, percentil):
    # Primero se crea una copia del DataFrame original y se representa el strip plot con seaborn
    df_plot = df.copy()
    
    # Generamos los subplots
    fig, axex = plt.subplots(1, 2, figsize = (10, 8), sharey = True)

    # Se establecen los colores para los boxplots y los puntos de cada método, y se configura el estilo de la gráfica con seaborn
    sns.set_theme(style = "ticks")
    box_colors = {"DESeq2": "#9ecae1", "limma": "#a1d99b", "edgeR": "#fdae6b"}
    point_colors = {"DESeq2": "#3182bd", "limma": "#31a354", "edgeR": "#e6550d"}
    
    # Se obtiene la lista de asignaciones únicas para iterar sobre ellas
    asignacion = df_plot["assignment"].unique().tolist()

    # Y similar a la representación anterior se representa esta vez el percentil
    for i, asing in enumerate(asignacion):
        ax = axex[i]
        # Pintamos el boxplot de fondo para cada método, con un color específico para cada uno
        sns.boxplot(data = df_plot[df_plot["assignment"] == asing], x = "db", y = "abs(FoldChange)_Percentil_sig", hue = "metodo", palette = box_colors, width = 0.6, dodge = True, fliersize = 0, boxprops = dict(alpha = 0.6), ax = ax)

        # Y pintamos los puntos por encima, con un color específico para cada método
        sns.stripplot(data = df_plot[df_plot["assignment"] == asing], x = "db", y = "abs(FoldChange)_Percentil_sig", hue = "metodo", palette = point_colors, dodge = True, jitter = 0.15, size = 5, alpha = 0.8, marker = "o", ax = ax)

        # Se añade el titulo a cada subplot y se eliminan marcas individuales de los ejes
        ax.set_title(f"{asing}", fontweight = "bold", pad = 12)
        ax.legend_.remove()
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Nombramos cada eje y la gráfica, añadiendo también una leyenda para los métodos
    fig.supxlabel("Base de Datos")

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[:3], labels[:3], loc = "upper left", ncol = 1, title = "Marco de Expresion Diferencial", title_fontsize = "11", fontsize = "9",  bbox_to_anchor = (1.02, 1), frameon = True)
    fig.suptitle("Efecto de la base de datos y normalizacion sobre el Fold Change", fontweight = "bold")
    fig.supylabel(f"Percentil {percentil} del abs(log2FoldChange) de miRNAs significativos", fontsize = "9")
    sns.despine()
    plt.tight_layout()
    plt.savefig("percentil_log2FoldChange.png", dpi = 300, bbox_inches = "tight")
    plt.close()
    return

def main():
    files, output_path = select_archives()
    if len(files) == 0:
        print("No se seleccionaron archivos. Saliendo del programa.")
        sys.exit(1)
    # Cambiamos el directorio de trabajo a la ruta de salida proporcionada por el usuario
    os.makedirs(output_path, exist_ok = True)
    os.chdir(output_path)
    percentil = 80
    df = data_frame(files, percentil)
    # Se crean los strip plots 
    strip_plot_padj(df, "miRNAs_significativos_DE")
    strip_plot_padj(df, "porcentaje_sig_estudio")
    strip_plot_padj(df, "porcentaje_sig_db")

    strip_plot_percentil_log2fc(df, percentil)


if __name__ == "__main__":
    main()