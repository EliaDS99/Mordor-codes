#!/bin/bash
#SBATCH --account IscrC_GONDOR
#SBATCH -p g100_usr_prod
#SBATCH --time 04:00:00
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --job-name=mordor_for_loop

# -----------------------------------------------
# Esempio di run_mordor.sh, basato su run_all.sh:
# accetta in input un file di testo (es. gal_pot_arepo.txt)
# contenente — riga per riga — tutti i nomi “gal_potupdated_XXXXX.hdf5”.
# Per ciascun nome, invoca: python3 mordor.py <galfile> --mode iso_sim --ShowPlots
# -----------------------------------------------

module load python/3.11.7--gcc--10.2.0

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR

if [ -z "$1" ]; then
  echo "Usage: sbatch run_mordor.sh nome_file_di_lista.txt"
  exit 1
fi

LISTFILE="$1"

# Controlla che il file LISTFILE esista
if [ ! -f "$LISTFILE" ]; then
  echo "ERROR: Il file '$LISTFILE' non esiste in $(pwd)"
  exit 1
fi

echo "DEBUG: Leggerò la lista di file HDF5 da: $LISTFILE"
echo "DEBUG: Contenuto di $LISTFILE:"
cat "$LISTFILE"
echo "DEBUG: Numero di righe: $(wc -l < "$LISTFILE")"
echo

# (Opzionale) Creo una cartella per salvare i log di Mordor
LOGDIR="logs_mordor"
mkdir -p "$LOGDIR"

# Loop su ogni riga del file di lista
while IFS= read -r galfile || [ -n "$galfile" ]; do
  # Salto righe vuote o commentate
  if [[ -z "$galfile" || "$galfile" =~ ^\# ]]; then
    echo "DEBUG: Riga vuota o commentata, salto."
    continue
  fi

  echo "-----------------------------------------------"
  echo "Processing: $galfile"
  echo "-----------------------------------------------"

  # Controlla che il file HDF5 esista
  if [ ! -f "$galfile" ]; then
    echo "ERROR: File HDF5 '$galfile' non trovato, salto."
    continue
  fi

  # Estrazione del nome base (senza estensione) per eventualmente rinominare output
  base=$(basename "$galfile" .hdf5)
  # Esempio: galfile="gal_potupdated_000004.hdf5" → base="gal_potupdated_000004"

  # Esegui Mordor con 4 CPU (o il numero che hai allocato)
  # Reindirizzo stdout e stderr in un log specifico per questa galassia
  python3 mordor2.py "$galfile" --mode iso_sim --ShowPlots \
    &> "${LOGDIR}/${base}.log"

  RET=$?
  if [ $RET -ne 0 ]; then
    echo "ERROR: mordor.py è terminato con codice $RET per $galfile"
    echo "  -> Controlla il log: ${LOGDIR}/${base}.log"
  else
    echo "OK: terminato $galfile, vedi log in ${LOGDIR}/${base}.log"
  fi

  # (Opzionale) Sposta i file .png generati in una cartella dedicata a questa galassia
  OUTDIR="results_${base}"
  mkdir -p "$OUTDIR"
  mv "${base}"_*.png "$OUTDIR" 2>/dev/null || true
  mv "${base}"_rbins.txt "$OUTDIR" 2>/dev/null || true
  echo
done < "$LISTFILE"

echo "DEBUG: Loop concluso su $LISTFILE"

