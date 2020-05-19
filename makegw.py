import pylab
from pycbc.waveform import get_td_waveform
import matplotlib.pyplot as plt
import numpy as np

import h5py


apx = 'SEOBNRv4'

RATE = 8192
dt = 1. / RATE

m_max = 75.
m_min = 5.
dm = 0.5

m_num = int((m_max-m_min)/dm+1)

################
GM0 = 6.674e-11 * 1.9891e30
c   = 299792458.
C1=5**(3/8.)/(8*np.pi) * ((c**3)/GM0)**(5./8) 

def time2lowfreq(m1,m2,t):
    Mcf= (m1*m2)**(-0.375) * (m1+m2)**(0.125)
    freq = C1 * Mcf * (t*1.4)**(-0.375)    ## here I add factor 1.4 just to get a longer waveform
    
    return freq
################

mA = []
mB = []



H5_FILE = 'bbh_%d_dm%1d.h5' %(RATE, dm)
#try:
f = h5py.File(H5_FILE, 'w', libver='latest')
#except:
#    f.close()
#    f = h5py.File(H5_FILE, 'w', libver='latest')
    
main_grep = f.create_group('/waveform')
main_grep.attrs['srate'] = RATE
main_grep.attrs['model'] = apx
main_grep.attrs['desc'] = 'Spinless BBH waveform model'
    

step = 0

for i in range(0, m_num):
    for j in range(i, m_num):
        ma = m_max - i*dm
        mb = m_max - j*dm
        
        f_low = time2lowfreq(ma, mb, 2.0)
        
        mA.append(ma)
        mB.append(mb)
        
        hp, hc = get_td_waveform(approximant=apx,
                    mass1=ma, mass2=mb, spin1z=0, delta_t=dt, f_lower=f_low)
        hp_len = np.array(hp).shape[0] - RATE -1
        hp_ = hp[hp_len:-1]
        hc_len = np.array(hc).shape[0] - RATE -1
        hc_ = hc[hc_len:-1]
        
        print('#%s with m1: %.1f, m2: %.1f, f_low = %.1f, t=[%f %f], len= %d'
              %(step,ma,mb,f_low,hp_.sample_times[0],hp_.sample_times[-1],np.array(hp_).shape[0]))

        gname = '%s' %(step)
        
        grp = main_grep.create_group(gname)
        grp.attrs['m'] = [ma, mb]
        grp.attrs['sz'] = [0,0]
        grp.attrs['F_low'] = f_low
        grp.create_dataset('hp', data=hp_, dtype='f')
        grp.create_dataset('hc', data=hc_, dtype='f')

        step = step + 1

f.close()

mA_ = np.array(mA)
mB_ = np.array(mB)
#print(np.array(mA).shape)
#print(np.array(mB).shape)

plot = False
if plot:
    plt.figure(figsize=(8,8))
    plt.plot(mA_, mB_, 'r.', markersize=2)
    plt.savefig('mass.pdf')
    plt.show()
        
        
plot_ = False
if plot_:
    plt.figure()
    plt.plot(hp_)
        
        
        
        
        
        
        
        
        
        
        
        
