'''
Make the plot of Lflare/Lkp vs Lhalpha/Lbol

Do we see coherence?
'''


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

flarefile = '../data/kic_lflare.csv'
fdata = pd.read_csv(flarefile)

'''
match Amy's M dwarf tables to the RA,Dec
Marcel will then use these to look for serendipitous X-ray detections
'''

Mfile = '../data/mcquillan/all_m.csv'
kic_m, kic_mper = np.loadtxt(Mfile, delimiter=',', usecols=(0,1), skiprows=1, unpack=True)

ra_m = np.zeros_like(kic_m)-99
de_m = np.zeros_like(kic_m)-99
gi_m = np.zeros_like(kic_m)-99
fl_m = np.zeros_like(kic_m)-99 # the Lfl/Lkp measurement!

for k in range(len(kic_m)):
    mtc = np.where((kic_m[k] == ddata['kicnum'].values))
    ra_m[k] = fdata['ra'].values[mtc][0]
    de_m[k] = fdata['dec'].values[mtc][0]
    fl_m[k] = fdata['LflLkep'].values[mtc][0]
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
hfile = '../data/lamost_match/amy_lamost_xmatch_extra.csv' # file w/ LAMOST and Hammer results
hdata = pd.read_csv(hfile)

ewha_m = np.zeros_like(kic_m)-99
ewhaerr_m = np.zeros_like(kic_m)-99
file_L = np.empty(len(kic_m), dtype='S30')

for k in range(len(kic_m)):
    mtc = np.where(ra_m[k] == hdata['col1'].values)
    if len(mtc[0])>0:
        ewha_m[k] = hdata['ewHa '].values[mtc][0]
        ewhaerr_m[k] = hdata['ewerr '].values[mtc][0]
        file_L[k] = hdata['eyename'].values[mtc][0]

# only the M dwarfs in our sample with valid EWHA measurements
ok = np.where((ewha_m > -99) & (ewhaerr_m < 100))
ok1 = np.where((ewha_m > -99) & (ewhaerr_m < 100) & (kic_mper > 0))
ok0 = np.where((ewha_m > -99) & (ewhaerr_m < 100) & (kic_mper <= 0))

# dfout = pd.DataFrame(data={'kicnum':kic_m[ok], 'Lfile':file_L[ok]})
# dfout.to_csv('lamost_matches.csv')


# http://adsabs.harvard.edu/abs/2014ApJ...795..161D
# file from here: https://figshare.com/articles/Chi_values_for_K_and_M_dwarfs_Table_8_in_Douglas_14_/1275226
chifile = '../data/chi_douglas2014.tsv'
gr, ri, logchi = np.loadtxt(chifile, unpack=True, usecols=(2,3,10))

gi_chi = gr + ri

# get the chi value for each star
logchi_m = np.interp(gi_m, gi_chi, logchi)

lha_lbol_m = ewha_m * logchi_m * 10e-5
lha_lbolerr_m = ewhaerr_m * logchi_m * 10e-5
######### /LAMOST

######### Leslie & Tessa
tfile = '../data/tessa_red_ew1.txt'

kic_t, ewha_t = np.loadtxt(tfile, usecols=(0,1), unpack=True, skiprows=1)

ewha_mt = np.zeros_like(kic_m)-99
for k in range(len(kic_t)):
    mtc = np.where((kic_t[k] == kic_m))
    if len(mtc[0])>0:
        ewha_mt[mtc] = ewha_t[k]

lha_lbol_mt = ewha_mt * logchi_m * 10e-5

okt = np.where((ewha_mt > -99))
okt1 = np.where((ewha_mt > -99) & (kic_mper > 0))
okt0 = np.where((ewha_mt > -99) & (kic_mper <= 0))

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
# plt.figure()
# plt.scatter(ewha_m[ok], np.log10(fl_m[ok]))
# plt.xlim(-3,5)
# plt.ylim(-9,-2)
# plt.xlabel('EW Halpha')
# plt.ylabel('log L$_{flare}$/L$_{kp}$')
# plt.savefig('figures/flare_vs_ewha.png')
# plt.close()

print('# LAMOST stars with Prot meas: ' + str(len(ok1[0])))
print('# LAMOST stars without Prot meas: ' + str(len(ok0[0])))

print('# Tessa stars with Prot meas: ' + str(len(okt1[0])))
print('# Tessa stars without Prot meas: ' + str(len(okt0[0])))

mm = 0
for k in range(len(ok)):
    xm = np.where((okt[0] == ok[0][k]))
    if len(xm[0])>0:
        mm = mm + 1
print('# of stars from LAMOST with Tessa data: '+str(mm))


plt.figure()
# LAMOST stars
plt.scatter(lha_lbol_m[ok], np.log10(fl_m[ok]), c='k')
plt.errorbar(lha_lbol_m[ok], np.log10(fl_m[ok]), xerr=lha_lbolerr_m[ok],
             fmt='o', color='k', marker='.')
# plt.scatter(lha_lbol_m[ok1], np.log10(fl_m[ok1]), c='k', marker='*', s=40)

# data for Tessa's stars
plt.scatter(lha_lbol_mt[okt], np.log10(fl_m[okt]), c='r', alpha=0.7)
# plt.scatter(lha_lbol_mt[okt1], np.log10(fl_m[okt1]), c='r', marker='*', s=40)

'''
# data for GJ 1243, GJ 1245A, GJ 1245B
plt.scatter([10**-3.56, 10**-4.14, 10**-3.97],
            [-3.78, -3.93, -4.00], marker='^', alpha=0.5)
'''

plt.xlim(-0.001, 0.003)
plt.ylim(-9,-2)
plt.xlabel(r'LH$\alpha$/L$_{bol}$')
plt.ylabel('log L$_{flare}$/L$_{kp}$')
plt.savefig('figures/flare_vs_lhalbol.png',dpi=150)
plt.close()


# plt.figure()
# plt.scatter(kic_mper[ok], lha_lbol_m[ok], c='k')
# plt.scatter(kic_mper[okt], lha_lbol_mt[okt], c='r')
# plt.show()
