#!/usr/bin/env Python

import pandas as pd
import re
import argparse
import os

# --- Configurar los argumentos ---
parser = argparse.ArgumentParser(description="Convertir FASTA de UniProt a TSV")
parser.add_argument("input", help="Archivo FASTA de entrada")
args = parser.parse_args()

# --- Generar nombre de salida automático ---
base, _ = os.path.splitext(args.input)
output_file = base + ".tsv"

# --- Expresión regular para parsear los encabezados ---
patron = re.compile(
    r"^>(?P<db>\w+)\|(?P<acceso>\w+)\|(?P<entry>\S+)\s+(?P<desc>.+?)\s+OS=(?P<organismo>.+?)\s+OX=(?P<ox>\d+)\s+GN=(?P<gen>\S+)\s+PE=(?P<pe>\d+)\s+SV=(?P<sv>\d+)"
)

data = []

# --- Leer el archivo de entrada ---
with open(args.input, "r", encoding="utf-8") as f:
    for linea in f:
        if linea.startswith(">"):
            m = patron.match(linea.strip())
            if m:
                data.append(m.groupdict())

# --- Guardar en TSV ---
df = pd.DataFrame(data)
df.to_csv(output_file, sep="\t", index=False)

print(f"✅ Archivo TSV generado automáticamente: {output_file}")

