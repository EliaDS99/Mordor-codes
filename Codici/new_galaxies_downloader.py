import os
import requests
import h5py
import numpy as np
import time

# Parameters
SimName = 'TNG100-1'
SnapNum = 99

ncut = 1e4  # Minimum number of star particles

# Stellar mass limits (in units of 1e10 Msun/h)
mass_min = 0.6774          # e.g. 1e10 Msun/h
mass_max = None           # Set to None to disable upper mass limit

# Maximum number of galaxies to download; set to None for no limit
max_galaxies = None       # e.g. 100 for the first 100 galaxies, or None for all

# Flags for each particle type; set include_potential to False to remove 'Potential'
include_potential = False
flags = {
    'gas': ['Coordinates', 'Masses', 'ParticleIDs', 'Velocities', 'Potential'],
    'dm': ['Coordinates', 'ParticleIDs', 'Velocities', 'Potential'],
    'stars': ['Coordinates', 'Masses', 'ParticleIDs', 'Velocities', 'Potential'],
    'bhs': ['Coordinates', 'Masses', 'ParticleIDs', 'Velocities', 'Potential']
}

# Remove 'Potential' if requested
if not include_potential:
    for key in flags:
        flags[key] = [f for f in flags[key] if f != 'Potential']

# Build params dictionary for requests
cutout_params = {ptype: ','.join(fields) for ptype, fields in flags.items()}

# Base API settings
baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key": "1701eec7f65350338e690222837ebc97"}  # Inserisci la tua vera API key

# Discover simulation and snapshot URLs
r = requests.get(baseUrl, headers=headers).json()
names = [sim['name'] for sim in r['simulations']]
i = names.index(SimName)
sim = requests.get(r['simulations'][i]['url'], headers=headers).json()
snaps = requests.get(sim['snapshots'], headers=headers).json()
snap = requests.get(snaps[SnapNum]['url'], headers=headers).json()

# Query subhalos with M_star > mass_min
search_query = f"?mass_stars__gt={mass_min}&filterFlag=True"
subs = snap['subhalos'] + search_query
subhalos = requests.get(subs, headers=headers).json()
ntot = subhalos['count']
print(f"Found {ntot} subhalos with M_star > {mass_min}")

# Retrieve all matching subhalo IDs
subhalos = requests.get(subs + f"&limit={ntot}", headers=headers).json()
ids = [entry['id'] for entry in subhalos['results']]

# Apply maximum galaxies limit if specified
if max_galaxies is not None:
    ids = ids[:max_galaxies]
    print(f"Limiting to first {len(ids)} subhalos")
else:
    print(f"No limit on number of subhalos; processing all {len(ids)}")

# Estimate total download size via GET with stream to get Content-Length
#print("Estimating total download size (fields: ", cutout_params, ")...")
#total_bytes = 0
#for idx, sub_id in enumerate(ids, start=1):
    # Every 500 entries, print a status message
#    if idx % 500 == 0 or idx == 1:
#        print(f"Sto contando: {idx}/{len(ids)} subhalos")
#    try:
#       url = snap['subhalos'] + str(sub_id) + "/cutout.hdf5"
#      r_stream = requests.get(url,
#                                params=cutout_params, headers=headers, stream=True, timeout=60)
#        r_stream.raise_for_status()
#        size = int(r_stream.headers.get('Content-Length', 0))
#        total_bytes += size
#        # Print size of this subhalo's cutout
#        print(f"  Subhalo {sub_id}: {size/1024/1024:.2f} MB")
#        r_stream.close()
#    except Exception:
#        continue
# After loop, print final count
#print(f"Sto contando: {len(ids)}/{len(ids)} subhalos")
#total_mb = total_bytes / (1024**2)
#total_gb = total_mb / 1024
#print(f"Total size to download: {total_gb:.2f} GB")

# Ask for confirmation
#start = input("Start download? (y/n): ")
#if start.lower() != 'y':
#    print("Download canceled.")
#    exit()

# Prepare download
downloaded_files = []
spinner = ['|', '/', '-', '\\']
spin_idx = 0

# Download loop: iterate over selected subhalo IDs
for sub_id in ids:
    # Spinner indicator
    print(f"\r{spinner[spin_idx]} Downloading galaxy {sub_id}", end='', flush=True)
    spin_idx = (spin_idx + 1) % len(spinner)
    try:
        url = snap['subhalos'] + str(sub_id)
        gal_info = requests.get(url, headers=headers).json()
        m_star = gal_info['mass_stars']
        # Check upper mass limit if specified
        if mass_max is not None and m_star > mass_max:
            continue
        if gal_info['len_stars'] <= ncut:
            continue
        # Download HDF5 cutout
        r_file = requests.get(url + "/cutout.hdf5",
                              params=cutout_params, headers=headers)
        r_file.raise_for_status()
        filename = r_file.headers['content-disposition'].split('filename=')[1]
        with open(filename, 'wb') as f:
            f.write(r_file.content)
        new_name = f"gal_{str(sub_id).zfill(6)}.hdf5"
        os.rename(filename, new_name)
        downloaded_files.append(new_name)
        # Update HDF5 header attributes for AREPO compatibility
        with h5py.File(new_name, 'r+') as f:
            header = f['Header'].attrs
            if 'NumPart_ThisFile' in header:
                header['NumPart_Total'] = header['NumPart_ThisFile']
                header['NumPart_Total_HighWord'] = np.zeros(6, dtype='i4')
                header['NumPart_ThisFile_HighWord'] = np.zeros(6, dtype='i4')
            for flag in ['Flag_Sfr', 'Flag_Cooling', 'Flag_StellarAge', 'Flag_Metals', 'Flag_Feedback', 'Flag_DoublePrecision']:
                header[flag] = np.int32(0)
    except Exception:
        print("\nOps... something gone wrong")
        break

# Write list of downloaded galaxy files to disk
if downloaded_files:
    ids_nums = [int(fn.split('_')[1].split('.')[0]) for fn in downloaded_files]
    first_id = min(ids_nums)
    last_id = max(ids_nums)
    output_list_file = f"logfile_{first_id}_{last_id}.txt"
else:
    output_list_file = "logfile_empty.txt"

with open(output_list_file, 'w') as f:
    for name in downloaded_files:
        f.write(name + "\n")
print(f"\nWrote list of downloaded galaxies to {output_list_file}")
