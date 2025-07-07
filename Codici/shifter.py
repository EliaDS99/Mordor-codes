import h5py
import numpy as np
import sys

if __name__ == "__main__":
	logfilename = str(sys.argv[1])

with open(logfilename, 'r') as logfile:
	filenames = logfile.read().splitlines()

for filename in filenames:

	f = h5py.File(filename, 'r+')

	L = f['/Header'].attrs['BoxSize']
	s = L/2. -np.array(f['/PartType4/Coordinates'][0])
	c = 'Coordinates'

	for g in f.keys():
		if g=='Header':
			continue

		if c in f[g].keys():
			for i in np.arange(3):
				f[g][c][:,i] += s[i]
				f[g][c][:,i] = np.where(f[g][c][:,i]<0, f[g][c][:,i]+L, f[g][c][:,i])
				f[g][c][:,i] = np.where(f[g][c][:,i]>L, f[g][c][:,i]-L, f[g][c][:,i])

	f.close()