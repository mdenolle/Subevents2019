## Packages required :

import numpy as np
import matplotlib.pyplot as plt
import sys,os, glob
import os.path
from scipy import stats
import matplotlib
from matplotlib import ticker
pi=np.pi


def gaussienne(scale,center,standard,time_vector) :
    '''
    Returns a Gaussian Function, centered on [center] and with a standard deviantion (related to width of the curve) of [std].
    
    Time_vector must be an array
    
    For more :
    
    - Width at mi-max is 2.355*[std].
     '''
    # print('Width at mi-max = ',2.355*standard)
    return (np.exp(-(time_vector-center)**2/(2*standard*standard))/(standard*np.sqrt(2*pi)))*scale


dt=0.0703 #  s,  sampling space of SCARDEC

path='./SIM_STF/*.dat'
quakelist=glob.glob(path)               # make a list of filenames


# variable initization
Ms = np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Tsub=np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Dsub=np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Mw = np.zeros(len(quakelist))
M0 = np.zeros(len(quakelist))
mo1 = np.zeros(len(quakelist))
MM0 = np.zeros(len(quakelist))
Td = np.zeros(len(quakelist))
err = np.zeros(len(quakelist))
Mw = np.zeros(len(quakelist))
Nsub = np.zeros(len(quakelist),dtype=np.int)



for i,filename in enumerate(quakelist):
    print(filename)
    piplot=np.loadtxt(filename,skiprows=1)
    opened=open(filename)
    opened.readline()
    # read STF:
    time=np.zeros(len(piplot));rate=np.zeros(len(piplot))            # initialize time and moment-rate vectors
    for ii in range(len(piplot)) : # read each time stamp
        time[ii]=(piplot[ii][0])
        rate[ii]=(piplot[ii][1]) 
    # find index of positive time and reliable amplitudes.
    I=np.where( (rate>=0.001*np.max(rate)) & (time>=0)) [0]
    Td[i]=(time[I[-1]])           # duration of quake
    M0[i]=(np.trapz(rate[I],x=time[I])) # moment calculated from integrating the STF
    rate=rate/M0[i]     # We normalize the STF
    rate0=rate # STF that will not undergo the Gaussian substractions

    sub=0  # initially, no peak detected => 0 peaks
    gauss_final=np.zeros(len(rate)) # Final  Gaussian-built STF, of the same size as 'rate'
    for el in I :                       # go through time
        if rate[el-1]<rate[el] and rate[el]>rate[el+1] and rate[el]>(0.1*max(rate0)) and time[el]>0  : # peak detection / default = 0.10 for the min. value of peak
            error0=1e99  # initial error for the grid fit
            std0=0.
            for std in np.linspace(0.01/2.335,300/2.335,700) : # grid fit
                gauss=gaussienne(rate[el],time[el],std,time)
                gauss=gauss*rate[el]/max(gauss) 
                error=np.sum((rate[el-5:el+5]-gauss[el-5:el+5])**2)
                if error<error0 :
                    std0=std
                    error0=error
                    gauss0=gauss
#                 # Computation of the event's magnitude given the subevent's magnitude
        # if duration is greater than 1s and shorter than entire source duration
            if std0>1/4 and 4*std0< 1.2*Td[i]:
                    gauss_final=gauss_final+gauss0 # sum up the subevent gaussian
                    rate=rate-gauss0    # make residual
                    sub+=1              # increment subevent

                    Ms[i,sub-1]=np.trapz(gauss0,x=time)*M0[i]  # store moment of each subevent
                    Tsub[i,sub-1]=time[el] # store time at which it occurs
                    Dsub[i,sub-1]=std0*2*np.sqrt(2*np.log(10)) # duration of subevent
     
    MM0[i]=(np.trapz(gauss_final,x=time)*M0[i] )    # recover reconstructed moment 
    Nsub[i]=(sub)                           # store number of subevents for that quake
    err[i]=(MM0[i]/M0[i])                  # store error between reconstructed and true moment.
    Mw[i]=(2/3*np.log10(M0[i])-6.07)        # store moment magnitude
    print(Mw[i])
# store in variable. 
np.savez('allvar_simulations',M0=M0,Ms=Ms,MM0=MM0,Nsub=Nsub,Tsub=Tsub,Dsub=Dsub,Td=Td,err=err,Mw=Mw)


data=np.load('allvar_simulations.npz')
M0=data['M0'];Ms=data['Ms']
MM0=data['MM0'];Nsub=data['Nsub'];Nsub=Nsub.astype(np.int)
Tsub=data['Tsub'];Dsub=data['Dsub']
Td=data['Td'];#FM=data['FM']
# depth=data['depth']
err=data['err']
Mw=data['Mw']


plt.semilogx(M0,err,'o')
plt.show()

# #2: scardec
mmm0=[];mmms=[];mmsub=[];fmm=[];nsuub=[]
for i in range(len(Nsub)):
    for ii in range(np.minimum(Nsub[i],len(Ms[i,:]))):
        mmm0.append(M0[i])
        mmms.append(Ms[i,ii])
        nsuub.append(Nsub[i])

##### FIGURE S5  ###############
plt.figure(figsize=(8,11))
II=np.where(Nsub!=0)[0]
print(len(II),len(quakelist))
Mwbin=np.linspace(6,9,4).astype(np.int)
M0bin=10**(3/2*np.linspace(6,9,4)+9.1)
ii=np.where(np.asarray(nsuub)>0)[0].astype(np.int)
plt.subplot(211)
plt.scatter(np.asarray(mmm0)[ii],np.asarray(nsuub)[ii],edgecolor=None,c='b',alpha=0.3,s=30)
plt.xscale('log')
plt.grid(True,linewidth=0.25)
plt.rcParams['axes.axisbelow'] = True
plt.xlabel('$M_0$ (Nm)',fontsize=14);plt.ylabel('Number of subevents',fontsize=14)
plt.ylim(0,30);#plt.xlim(1E17,1E23)
plt.title('a)',loc='left')
ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)


##### FIGURE 5b: PLOT MOMENT2MOMENT ###############
M0bin1=10**(3/2*np.linspace(6,9,9)+9.1)
# find the median
medMs=[];medM0=[]
for i in range(len(M0bin1)):
        ik=np.where((mmm0<=M0bin1[i]) & (mmm0>M0bin1[i-1]))[0]
        medMs.append(10**(np.median(np.log10(np.asarray(mmms)[ik]))))
        medM0.append(10**(np.mean(np.log10(M0bin1[i-1:i]))))


plt.subplot(212)
# SCARDEC
plt.grid(True,linewidth=0.25)
plt.scatter(mmm0,mmms,c='b',edgecolor='k',alpha=0.3,s=2)
plt.xscale('log');plt.yscale('log')
plt.xlabel('$M_0$ (Nm)',fontsize=14)
plt.ylabel('$M_S$ (Nm)',fontsize=14)
plt.title('b) ',loc='left',fontsize=14)

ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)

plt.show()