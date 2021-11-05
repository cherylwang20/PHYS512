import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
#plt.ion()

def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    #gpsStart=meta['GPSstart'].value
    gpsStart=meta['GPSstart'][()]
    #print meta.keys()
    #utc=meta['UTCstart'].value
    utc=meta['UTCstart'][()]
    #duration=meta['Duration'].value
    duration=meta['Duration'][()]
    #strain=dataFile['strain']['Strain'].value
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc



#fnames=glob.glob("[HL]-*.hdf5")
#fname=fnames[0]
H=['H-H1_LOSC_4_V2-1126259446-32.hdf5','H-H1_LOSC_4_V2-1135136334-32.hdf5','H-H1_LOSC_4_V1-1167559920-32.hdf5','H-H1_LOSC_4_V2-1128678884-32.hdf5']
L = ['L-L1_LOSC_4_V2-1126259446-32.hdf5','L-L1_LOSC_4_V2-1135136334-32.hdf5','L-L1_LOSC_4_V1-1167559920-32.hdf5','L-L1_LOSC_4_V2-1128678884-32.hdf5']
strain_H = [0]*len(H)
strain_L = [0]*len(H)
utc_H = [0]*len(H)
utc_L = [0]*len(H)
dt_H = [0]*len(H)
dt_L = [0]*len(H)
for i in range(len(H)):
    print('reading file ',H[i],'and',L[i])
    strain_H[i],dt_H,utc_H[i]= read_file(H[i])
    strain_L[i], dt_L, utc_L[i] = read_file(L[i])

#th,tl=read_template('GW150914_4_template.hdf5')
GW=['GW150914_4_template.hdf5','GW151226_4_template.hdf5','GW170104_4_template.hdf5','LVT151012_4_template.hdf5']



th = [0]*len(GW)
tl = [0]*len(GW)
for i in range(len(GW)):
    th[i],tl[i]=read_template(GW[i])


#spec,nu=measure_ps(strain,do_win=True,dt=dt,osamp=16)
#strain_white=noise_filter(strain,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)

#th_white=noise_filter(th,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)
#tl_white=noise_filter(tl,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)


#matched_filt_h=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(th_white)))
#matched_filt_l=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(tl_white)))




#copied from bash from class
# strain2=np.append(strain,np.flipud(strain[1:-1]))
# tobs=len(strain)*dt
# k_true=np.arange(len(myft))*dnu
