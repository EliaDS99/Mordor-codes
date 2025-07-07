import pandas as pd
import matplotlib.pyplot as plt

# 1) Carica il catalogo generato
df = pd.read_csv('morphology_catalog.txt', 
                 delim_whitespace=True, 
                 comment='#')

# 2) Calcola il rapporto D/T = (Mthin + Mthick + Mpbulge) / Mstar
df['D_T'] = (df['Mthin'] + df['Mthick'] + df['Mpbulge']) / df['Mstar']

# 3) Filtra eventuali NaN
d_t = df['D_T'].dropna()

# 4) Costruisci l’istogramma
plt.figure(figsize=(6,4))
plt.hist(d_t, bins=30, edgecolor='k')
plt.xlabel('Disk–to–Total mass ratio $D/T$')
plt.ylabel('Number of galaxies')
plt.title('Distribution of $D/T$ across all galaxies')
plt.tight_layout()

# 5) Salva la figura
plt.savefig('D_T_histogram.png')

