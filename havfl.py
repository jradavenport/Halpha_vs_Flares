'''
Make the plot of Lflare/Lkp vs Lhalpha/Lbol

Do we see coherence?
'''


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

flarefile = 'kic_lflare.csv'
radecfile = 'kic_radec.csv'
fdata = pd.read_csv(flarefile)
ddata = pd.read_csv(radecfile)


'''
match Amy's M dwarf tables to the RA,Dec
Marcel will then use these to look for serendipitous X-ray detections
'''

Mfile = 'amy/all_m.txt'
kic_m = np.loadtxt(Mfile, delimiter=',', usecols=(0,), skiprows=1, unpack=True)

ra_m = np.zeros_like(kic_m)-99
de_m = np.zeros_like(kic_m)-99
fl_m = np.zeros_like(kic_m)-99 # the Lfl/Lkp measurement!

for k in range(len(kic_m)):
    mtc = np.where((kic_m[k] == ddata['kicnum'].values))
    ra_m[k] = ddata['ra'].values[mtc][0]
    de_m[k] = ddata['dec'].values[mtc][0]
    fl_m[k] = fdata['LflLbol'].values[mtc][0]

### this was ouput, sent to Vizier's Xmatch tool to match to LAMOST
# dfout = pd.DataFrame(data={'ra':ra_m,'dec':de_m,'kicnum':kic_m})
# dfout.to_csv('allM_radec.csv')

### The object ID's from LAMOST that had matches in Amy's sample were then sent
### to the LAMOST DR1 site, and spectra retrieved (lamost_match/fitspng1938584811/)

### With LAMOST spectra in hand, I then ran Hammer to get Halpha

### now i need to re-match the Hammer outputs to Amy's KIC sample...
hfile = 'lamost_match/amy_lamost_xmatch_extra.csv' # file w/ LAMOST and Hammer results
hdata = pd.read_csv(hfile)

ewha_m = np.zeros_like(kic_m)-99
ewhaerr_m = np.zeros_like(kic_m)-99
for k in range(len(kic_m)):
    mtc = np.where(ra_m[k] == hdata['col1'].values)
    if len(mtc[0])>0:
        ewha_m[k] = hdata['ewHa '].values[mtc][0]
        ewhaerr_m[k] = hdata['ewerr '].values[mtc][0]

ok = np.where((ewha_m > -99))


plt.figure()
plt.scatter(ewha_m[ok], np.log10(fl_m[ok]))
plt.xlim(-3,5)
plt.ylim(-9,-2)
plt.xlabel('EW Halpha')
plt.ylabel('log L$_{flare}$/L$_{kp}$')
plt.savefig('flare_vs_ewha.png')
plt.close()




