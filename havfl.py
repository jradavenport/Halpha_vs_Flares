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
gi_m = np.zeros_like(kic_m)-99
fl_m = np.zeros_like(kic_m)-99 # the Lfl/Lkp measurement!

for k in range(len(kic_m)):
    mtc = np.where((kic_m[k] == ddata['kicnum'].values))
    ra_m[k] = ddata['ra'].values[mtc][0]
    de_m[k] = ddata['dec'].values[mtc][0]
    fl_m[k] = fdata['LflLbol'].values[mtc][0]
    gi_m[k] = fdata['giclr'].values[mtc][0]

### this was ouput, sent to Vizier's Xmatch tool (http://cdsxmatch.u-strasbg.fr/xmatch#tab=xmatch&)
### to match to LAMOST (V/146/dr1)
# dfout = pd.DataFrame(data={'ra':ra_m,'dec':de_m,'kicnum':kic_m})
# dfout.to_csv('allM_radec.csv')

### The object ID's from LAMOST that had matches in Amy's sample were then sent
### to the LAMOST DR1 site (http://dr1.lamost.org/q),
### and spectra retrieved (lamost_match/fitspng1938584811/)

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

# only the M dwarfs in our sample with valid EWHA measurements
ok = np.where((ewha_m > -99) & (ewhaerr_m < 100))

# http://adsabs.harvard.edu/abs/2014ApJ...795..161D
# file from here: https://figshare.com/articles/Chi_values_for_K_and_M_dwarfs_Table_8_in_Douglas_14_/1275226
chifile = 'chi_douglas2014.tsv'
gr, ri, logchi = np.loadtxt(chifile, unpack=True, usecols=(2,3,10))

gi_chi = gr + ri

# get the chi value for each star
logchi_m = np.interp(gi_m, gi_chi, logchi)

lha_lbol_m = ewha_m * logchi_m * 10e-5
lha_lbolerr_m = ewhaerr_m * logchi_m * 10e-5
######### /LAMOST

######### Leslie & Tessa
tfile = 'tessa_red_ew1.txt'

kic_t, ewha_t = np.loadtxt(tfile, usecols=(0,1), unpack=True, skiprows=1)

ewha_mt = np.zeros_like(kic_m)-99
for k in range(len(kic_t)):
    mtc = np.where((kic_t[k] == kic_m))
    if len(mtc[0])>0:
        ewha_mt[mtc] = ewha_t[k]

lha_lbol_mt = ewha_mt * logchi_m * 10e-5

okt = np.where((ewha_mt > -99))
######### /Leslie & Tessa

# Need to make sure: are the classic Kepler flare M dwarfs (GJ 1243, GJ 1245AB)
# in here? I don't think so yet... but we have good Halpha constraints!

'''
 GJ   Lfl/Lkp  Lha/Lbol
1243    −3.78   −3.56
1245A   −3.93   −4.14
1245B   −4.00   −3.97
'''





### Plots...!
plt.figure()
plt.scatter(ewha_m[ok], np.log10(fl_m[ok]))
plt.xlim(-3,5)
plt.ylim(-9,-2)
plt.xlabel('EW Halpha')
plt.ylabel('log L$_{flare}$/L$_{kp}$')
plt.savefig('flare_vs_ewha.png')
plt.close()


plt.figure()
plt.scatter(lha_lbol_m[ok], np.log10(fl_m[ok]), c='k')
plt.errorbar(lha_lbol_m[ok], np.log10(fl_m[ok]), xerr=lha_lbolerr_m[ok],
             fmt='o', color='k', marker='.')
plt.scatter(lha_lbol_mt[okt], np.log10(fl_m[okt]), c='r')
plt.scatter([10**-3.56, 10**-4.14, 10**-3.97],
            [-3.78, -3.93, -4.00], marker='^', s=40)
plt.xlim(-0.001, 0.003)
plt.ylim(-9,-2)
plt.xlabel(r'LH$\alpha$/L$_{bol}$')
plt.ylabel('log L$_{flare}$/L$_{kp}$')
plt.savefig('flare_vs_lhalbol.png',dpi=300)
plt.close()
