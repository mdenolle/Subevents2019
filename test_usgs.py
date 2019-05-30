## DISCLAIMER

# This program allows for the peak detection, gaussian fitting and parameters extractions of the subevents, as performed in Danre et al. (2019)
# Here we take the example of all Source Time Functions provided by Hayes (2017)
# The original version of the code was done by Philippe Danre.
# The code has been modified by Marine Denolle, last version June 1st 2019 (mdenolle@fas.harvard.edu)
# this is the exact same code as in subevent_GRL_script.py except that we do not look into the focal mechanism information and focus on the USGS database.
##
## import modules

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

path='./USGS_STFS/*.newstf'
quakelist=glob.glob(path)
M0=[];depth=[];FM=[];MM0=[];Nsub=[];dip=[];Td=[];mo1=[];err=[]
Mw=[]
# variable initization
Ms = np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Tsub=np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Dsub=np.zeros(shape=(len(quakelist),40),dtype=np.float64)
Mw = np.zeros(len(quakelist))
depth = np.zeros(len(quakelist))
dip = np.zeros(len(quakelist))
FM = np.zeros(len(quakelist))
M0 = np.zeros(len(quakelist))
mo1 = np.zeros(len(quakelist))
MM0 = np.zeros(len(quakelist))
Td = np.zeros(len(quakelist))
err = np.zeros(len(quakelist))
Mw = np.zeros(len(quakelist))
Nsub = np.zeros(len(quakelist),dtype=np.int)
for i,filename in enumerate(quakelist):
    
    print(i,filename)
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
    I = I[0:-1]
    Td[i]=(time[I[-1]])           # duration of quake
    M0[i]=(np.trapz(rate[I],x=time[I])) # moment calculated from integrating the STF
    rate=rate/M0[i]     # We normalize the STF
    rate0=rate # STF that will not undergo the Gaussian substractions

    sub=0  # initially, no peak detected => 0 peaks
    gauss_final=np.zeros(len(rate)) # Final  Gaussian-built STF, of the same size as 'rate'
    for el in I :                       # go through time
        if rate[el-1]<rate[el] and rate[el]>rate[el+1] and rate[el]>(0.1*max(rate0)) and time[el]>0  : # peak detection / default = 0.10 for the min. value of peak
            error0=1e99  # initial error for the grid fit
            std0=0.;gauss=0
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
     

    # one of the events is weird:
    if (2/3*np.log10(M0[i])-6.07 <= 7) and (sub>10):
        continue
    MM0[i]=(np.trapz(gauss_final,x=time)*M0[-1] )    # recover reconstructed moment 
    Nsub[i]=(sub)                           # store number of subevents for that quake
    err[i]=(MM0[i]/M0[i])                  # store error between reconstructed and true moment.
    Mw[i]=(2/3*np.log10(M0[i])-6.07)   
np.savez('allvar_USGS',M0=M0,Ms=Ms,MM0=MM0,Nsub=Nsub,Tsub=Tsub,Dsub=Dsub,Td=Td,FM=FM,depth=depth,err=err,Mw=Mw)



#exit()
data=np.load('allvar_USGS.npz')
M0=data['M0'];Ms=data['Ms']
MM0=data['MM0'];Nsub=data['Nsub']
Tsub=data['Tsub'];Dsub=data['Dsub']
Td=data['Td'];FM=data['FM']
depth=data['depth']
err=data['err']
Mw=data['Mw']


mmm0=[];mmms=[];mmsub=[];fmm=[];nsuub=[]
for i in range(len(quakelist)):
    for ii in range(Nsub[i]):
        mmm0.append(M0[i])
        mmms.append(Ms[i,ii])
        mmsub.append(Nsub[i])
        fmm.append(FM[i])
        nsuub.append(Nsub[i])



####### FIGURE S1 ##################
plt.figure(figsize=(8,11))
II=np.where(Nsub!=0)[0]
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(M0[II]),Nsub[II])
print(slope,intercept,r_value)
Mwbin=np.linspace(6,9,4).astype(np.int)
M0bin=10**((3/2)*(Mwbin+6.07))
ii=np.where(np.asarray(nsuub)>0)[0].astype(np.int)
plt.subplot(211)
plt.scatter(np.asarray(mmm0)[ii],np.asarray(mmsub)[ii],edgecolor=None,c='b',alpha=0.3,s=30)
plt.text(1E18,4,'slope='+str(round(slope*10)/10),color='r',fontsize=14,rotation=10)
plt.xscale('log')
plt.grid(True,linewidth=0.25)
plt.rcParams['axes.axisbelow'] = True
plt.xlabel('$M_0$ (Nm)',fontsize=14);plt.ylabel('Number of subevents',fontsize=14)
plt.ylim(0,40);plt.xlim(1E17,1E23)
plt.semilogx(M0,intercept + slope*np.log10(M0),linewidth=3,color='r')
plt.title('a)',loc='left')
ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)


Mwbin=np.linspace(6,9,4).astype(np.int)
M0bin=10**((3/2)*(Mwbin+6.07))
Mwbin1=np.linspace(5,9,18)
M0bin1=10**((3/2)*(Mwbin1+6.07))
# find the median
medMs=[];medM0=[]
for i in range(3,len(Mwbin1)):
        ik=np.where((mmm0<=M0bin1[i]) & (mmm0>M0bin1[i-1]))[0]
        medMs.append(10**(np.median(np.log10(np.asarray(mmms)[ik]))))
        medM0.append(10**(np.mean(np.log10(M0bin1[i-1:i]))))


plt.subplot(212)
plt.scatter(mmm0,mmms,c='b',edgecolor='k',alpha=0.3)
plt.xscale('log');plt.yscale('log')
plt.xlim(1E17,1E23);plt.ylim(1E15,1E23)
plt.xlabel('$M_0$ (Nm)',fontsize=14)
plt.ylabel('$M_S$ (Nm)',fontsize=14)
plt.grid(True,linewidth=0.25)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mmm0),np.log10(mmms))

print('regression for scaling with moments')
print(slope,intercept,r_value)
plt.text(1E21,1E16,'slope='+str(round(slope*100)/100),color='g',fontsize=14)
plt.loglog(mmm0,np.power(10,intercept + slope*np.log10(mmm0)),linewidth=3,color='g')
plt.title('b) ',loc='left',fontsize=14)
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23]),linewidth=1,color='r')
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23])/10,linewidth=1,color='r')
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23])/100,linewidth=1,color='r')

plt.scatter(medM0,medMs,marker='s',c='orange',edgecolor='k',s=40,zorder=10)

plt.text(1E18,7E18,'r = 1',color='r',fontsize=14,rotation=25)
plt.text(1E18,5E17,'r = 10',color='r',fontsize=14,rotation=25)
plt.text(1E18,7E16,'r = 100',color='r',fontsize=14,rotation=25)
ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)

plt.savefig('FigureS1.pdf')
plt.show()