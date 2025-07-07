#!/usr/bin/env bash
#SBATCH --job-name=plotter2
#SBATCH --output=output%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00

# Uso:
#   sbatch run_plotter2.sh galaxy_list.txt
#
# galaxy_list.txt deve contenere una riga per file, es:
#   gal_000000.hdf5
#   gal_000001.hdf5
#   ...

INPUT_LIST="$1"

# Parametri colormap e fattore di scala width = FACTOR * hmr
VMIN=1e2
VMAX=1e4
FACTOR=4

if [[ -z "$INPUT_LIST" ]]; then
  echo "Errore: specifica il file di input (lista di galassie)."
  echo "Uso: sbatch run_graf_tomm.sh galaxy_list.txt"
  exit 1
fi

echo "Avvio run_plotter2.sh su $INPUT_LIST"
while read -r galfile; do
  # estraggo l'indice dal nome, es: gal_000123.hdf5 → 000123
  index=$(basename "$galfile" .hdf5 | sed -E 's/.*_([0-9]+)$/\1/')
  echo -e "\n→ Lancio plotter2.py per index=$index (file=$galfile)"
  python3 plotter2.py "$index" --v "$VMIN" "$VMAX" "$FACTOR"
done < "$INPUT_LIST"

echo -e "\n→ Esecuzione completata su tutte le galassie."

