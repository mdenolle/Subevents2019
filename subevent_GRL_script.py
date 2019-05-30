## DISCLAIMER

# This program allows for the peak detection, gaussian fitting and parameters extractions of the subevents, as performed in Danre et al. (2019)
# Here we take the example of all Source Time Functions provided by SCARDEC (not sorted by focal mechanism, see the link : http://scardec.projects.sismo.ipgp.fr/). 

# For information about SCARDEC database, see Vallee et al. 2011 and Vallee et Douet 2016
# The original version of the code was done by Philippe Danre.
# The code has been modified by Marine Denolle, last version June 1st 2019 (mdenolle@fas.harvard.edu)
##
## import modules

import numpy as np
import matplotlib.pyplot as plt
import sys,os, glob
import os.path
import matplotlib.pyplot as plt
from scipy import stats
pi=np.pi




 # define a guaussian function
def gaussienne(scale,center,standard,time_vector) :
    '''
    Returns a Gaussian Function, centered on [center] and with a standard deviantion (related to width of the curve) of [std].
    
    Time_vector must be an array
    
    For more :
    
    - Width at half maximum is 2.355*[std]. Provides an useful time scale .
     '''
    return (np.exp(-(time_vector-center)**2/(2*standard*standard))/(standard*np.sqrt(2*pi)))*scale




dt=0.0703 #  s,  sampling space of SCARDEC


# read list of the ASCII files that have the SCARDEC database.
path='/Users/marinedenolle/Dropbox/SOURCE/SCARDEC/FCTS/*/*moy*'
quakelist=glob.glob(path)               # make a list of filenames


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

for i,filename in enumerate(quakelist): # loop through each quake

    # read SCARDEC files
    piplot=np.loadtxt(filename,skiprows=2)
    opened=open(filename)
    opened.readline()
    list=opened.readline().split(' ')   # read lines of parameters
    mo1[i]=(float(list[1]))          # extract moment
    depth[i]=(float(list[0]))        # extract depth
    dip[i]=(float(list[4]))          # extract dip
    r1=float(list[5]);r2=float(list[8]) # extract rake
     # use Shearer et al, 2006 to parameterize the focal mechanism type.
    if abs(r1)>90:r1=(180-abs(r1))*(r1/abs(r1))
    if abs(r2)>90:r2=(180-abs(r2))*(r2/abs(r2))
    if abs(r1)<abs(r2):
        FM[i]=r1/90
    else:
        FM[i]=r2/90

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
     
    MM0[i]=(np.trapz(gauss_final,x=time)*M0[-1] )    # recover reconstructed moment 
    Nsub[i]=(sub)                           # store number of subevents for that quake
    err[i]=(MM0[i]/M0[i])                  # store error between reconstructed and true moment.
    Mw[i]=(2/3*np.log10(M0[i])-6.07)        # store moment magnitude

# store in variable. 
np.savez('allvar',M0=M0,Ms=Ms,MM0=MM0,Nsub=Nsub,Tsub=Tsub,Dsub=Dsub,Td=Td,FM=FM,depth=depth,err=err,Mw=Mw)


# the following section plots figures in the main paper. you can stop here.
# for questions regarding the plotting, contact mdenolle@fas.harvard.edu (and check stackoverflow...)



## READ SCARDEC FIT
data=np.load('allvar.npz')
M0=data['M0'];Ms=data['Ms']
MM0=data['MM0'];Nsub=data['Nsub'];Nsub=Nsub.astype(np.int)
Tsub=data['Tsub'];Dsub=data['Dsub']
Td=data['Td'];FM=data['FM']
depth=data['depth']
err=data['err']
Mw=data['Mw']

## READ SIMULATION FIT (done with test_simulations.py)
data=np.load('allvar_simulations_G8.npz')
M0_sim=data['M0'];Ms_sim=data['Ms']
MM0_sim=data['MM0'];Nsub_sim=data['Nsub'];Nsub_sub=Nsub_sim.astype(np.int)
Tsub_sim=data['Tsub'];Dsub_sim=data['Dsub']
Td_sim=data['Td']

# indexes of those that did not rupture beyond the fault.
fidindez='/Users/marinedenolle/Dropbox/GROUP_PROJECTS/DANRE_DENOLLE_STRIKE_SLIP/good_indx_8.dat'
ID=open(fidindez,'r').readlines()
indx=np.zeros(len(ID),dtype=np.int)
for ii,i1 in enumerate(ID):
        indx[ii]=int(str(i1.split('\n')[0]))


# rearrange the variables to flatten with subevents:
# #1 simulations
mmm0_sim=[];mmms_sim=[];mmsub_sim=[]
for i in range(len(M0_sim)):
    if len(np.where(i==indx)[0])>=1: # only keep those that fit within the fault.
        for ii in range(Nsub_sim[i]):
            mmm0_sim.append(M0_sim[i])
            mmms_sim.append(Ms_sim[i,ii])
            mmsub_sim.append(Nsub_sim[i])
            
# #2: scardec
mmm0=[];mmms=[];mmsub=[];fmm=[];nsuub=[]
for i in range(len(quakelist)):
    for ii in range(Nsub[i]):
        mmm0.append(M0[i])
        mmms.append(Ms[i,ii])
        mmsub.append(Nsub[i])
        fmm.append(FM[i])
        nsuub.append(Nsub[i])
        


##### FIGURE 2a: PLOT SCARDEC SUBEVENTS ###############
plt.figure(figsize=(11,4.5))
II=np.where(Nsub!=0)[0]
II2=np.where( (FM[II]>=-0.5) & (FM[II]<=0.5)) [0]
II3=np.where( (FM[II]<-0.5) | (FM[II]>0.5)) [0]
II21=np.where( (FM[II]>=-0.5) & (FM[II]<=0.5) & (M0[II]>1.68E+19)) [0]
II22=np.where( (FM[II]>=-0.5) & (FM[II]<=0.5) & (M0[II]<=1.68E+19)) [0]
II32=np.where( ((FM[II]<-0.5) | (FM[II]>0.5)) & ( (M0[II]>=5E20)&(M0[II]<=2E22)) ) [0]
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(M0[II[II2]]),Nsub[II[II2]]) # all strike slip
slope0, intercept0, r_value, p_value, std_err = stats.linregress(np.log10(M0[II[II21]]),Nsub[II[II21]]) # all ss greater than 6.5
slope02, intercept02, r_value, p_value, std_err = stats.linregress(np.log10(M0[II[II22]]),Nsub[II[II22]]) #all ss smaller than 6.5
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.log10(M0[II[II3]]),Nsub[II[II3]]) # all dip slip
slope12, intercept12, r_value1, p_value1, std_err1 = stats.linregress(np.log10(M0[II[II32]]),Nsub[II[II32]]) # all greater than 7.8 but tohoku


Mwbin=np.linspace(6,9,4).astype(np.int)
M0bin=10**((3/2)*(Mwbin+6.07))
ii=np.where(np.asarray(nsuub)>0)[0].astype(np.int)
ii1=np.where(np.asarray(fmm)>=-0.5 )[0]
ii1=np.where( (np.asarray(fmm)[ii]>=-0.5) & (np.asarray(fmm)[ii]<=0.5)) [0]
ii1=ii[ii1]
ii2=np.where( (np.asarray(fmm)[ii]<-0.5) | (np.asarray(fmm)[ii]>0.5)) [0]
ii2=ii[ii2]
plt.subplot(121)
plt.scatter(np.asarray(mmm0)[ii2],np.asarray(mmsub)[ii2],edgecolor=None,c='b',alpha=0.3,s=30)
plt.scatter(np.asarray(mmm0)[ii1],np.asarray(mmsub)[ii1],edgecolor=None,c='r',alpha=0.1,s=24)
plt.text(3E21,15.5,'slope='+str(round(slope*10)/10),color='r',fontsize=14,rotation=25)
plt.text(3E21,5.5,'slope='+str(round(slope1*10)/10),color='b',fontsize=14,rotation=10)
plt.xscale('log')
plt.grid(True,linewidth=0.25)
plt.rcParams['axes.axisbelow'] = True
plt.xlabel('$M_0$ (Nm)',fontsize=14);plt.ylabel('Number of subevents',fontsize=14)
plt.ylim(0,30);plt.xlim(1E17,1E23)
plt.semilogx(M0,intercept + slope*np.log10(M0),linewidth=3,color='r')
plt.semilogx(M0,intercept1 + slope1*np.log10(M0),linewidth=3,color='b')
plt.title('a)',loc='left')
ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)



#### FIGURE 2b : PLOT SIMULATED SUBEVENTS ###############
plt.subplot(122)
plt.scatter(M0_sim,Nsub_sim,edgecolor=None,c='b',alpha=0.5,s=20)
plt.xscale('log')
plt.ylim((0,30))
plt.grid(True,linewidth=0.25)
plt.title('b)',loc='left')
plt.yticks([5,10,15,20,25,30])
plt.xlabel('$M_0$ (Nm/m)',fontsize=14);plt.ylabel('Number of subevents',fontsize=14)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mmm0_sim),mmsub_sim)
plt.semilogx(M0_sim,intercept + slope*np.log10(M0_sim),linewidth=3,color='g')
plt.text(2E14,6,'slope='+str(round(slope*10)/10),color='g',fontsize=14,rotation=10)
plt.savefig('Figure2.pdf')

##### FIGURE 3: PLOT MOMENT2MOMENT ###############
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


fig=plt.figure(figsize=(11,4.5))
plt.subplot(121)
# SCARDEC
plt.grid(True,linewidth=0.25)
# plt.axisbelow(True)
plt.scatter(mmm0,mmms,c='b',edgecolor='k',alpha=0.3)
plt.xscale('log');plt.yscale('log')
plt.xlim(1E17,1E23);plt.ylim(1E15,1E23)
plt.xlabel('$M_0$ (Nm)',fontsize=14)
plt.ylabel('$M_S$ (Nm)',fontsize=14)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mmm0),np.log10(mmms))

print('regression for scaling with moments')
print(slope,intercept,p_value)
plt.text(1E21,1E16,'slope='+str(round(slope*100)/100),color='g',fontsize=14)
plt.loglog(mmm0,np.power(10,intercept + slope*np.log10(mmm0)),linewidth=3,color='g')
plt.title('a) ',loc='left',fontsize=14)
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23]),linewidth=1,color='r')
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23])/10,linewidth=1,color='r')
plt.loglog(np.array([1E17,1E23]),np.array([1E17,1E23])/100,linewidth=1,color='r')

plt.scatter(medM0,medMs,marker='s',c='orange',edgecolor='k',s=40,zorder=10)

plt.text(1E22,1.5E22,'r = 1',color='r',fontsize=14,rotation=35)
plt.text(1E22,2E21,'r = 10',color='r',fontsize=14,rotation=35)
plt.text(1E22,4E20,'r = 100',color='r',fontsize=14,rotation=35)
ax2=plt.twiny()
ax2.set_xticks(np.array((np.log10(M0bin)-17)/6))
ax2.set_xticklabels(Mwbin.astype(str))
ax2.set_xlabel('$M_W$')
ax2.xaxis.set_label_coords(1, 1.05)



M0bin1=np.logspace(14.5,18,8)
# find the median
medMs=[];medM0=[]
for i in range(1,len(M0bin1)):
        ik=np.where((mmm0_sim<=M0bin1[i]) & (mmm0_sim>M0bin1[i-1]))[0]
        medMs.append(10**(np.median(np.log10(np.asarray(mmms_sim)[ik]))))
        medM0.append(10**(np.mean(np.log10(M0bin1[i-1:i]))))


plt.subplot(122)
plt.scatter(mmm0_sim,mmms_sim,edgecolor='k',c='b',alpha=0.5,s=20)
plt.xscale('log');plt.yscale('log')
plt.grid(True,linewidth=0.25)
plt.title('b)',loc='left')
plt.xlabel('$M_0$ (Nm/m)',fontsize=14);plt.ylabel('$M_S$ (Nm/m)',fontsize=14)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(mmm0_sim),np.log10(mmms_sim))
print('regression of moments for simulated')
print(slope,intercept)
plt.text(1E16,1E12,'slope='+str(round(slope*10)/10),color='g',fontsize=14)
plt.loglog(M0_sim,10**(intercept + slope*np.log10(M0_sim)),linewidth=3,color='g')
plt.loglog(np.array([1E14,1E17]),np.array([1E14,1E17]),linewidth=1,color='r')
plt.loglog(np.array([1E14,1E17]),np.array([1E14,1E17])/10,linewidth=1,color='r')
plt.loglog(np.array([1E14,1E17]),np.array([1E14,1E17])/100,linewidth=1,color='r')
plt.text(1E14,1E14,'r = 1',color='r',fontsize=14,rotation=30)
plt.text(1E14,1.5E13,'r = 10',color='r',fontsize=14,rotation=30)
plt.text(1E14,2E12,'r = 100',color='r',fontsize=14,rotation=30)
plt.scatter(medM0,medMs,marker='s',c='orange',edgecolor='k',s=40,zorder=10)
plt.savefig('Figure3.pdf')


##### FIGURE S4 ###############
fig=plt.figure(figsize=(8,10))
for i in range(len(quakelist)):
    if Nsub[i]==0:
        continue
    ibig=np.argmax(Ms[i,0:Nsub[i]])
    plt.subplot(311)
    plt.loglog(M0[i],Ms[i,0],'bo',markeredgecolor='black',markeredgewidth=0.2,alpha=0.5)
    plt.rcParams.update({'font.size': 14})
    plt.loglog(M0[i],Ms[i,ibig],'ro',markeredgecolor='black',markeredgewidth=0.2,alpha=0.5)
    plt.grid(True,linewidth=0.25)
    plt.title('a)',loc='left')
    plt.ylabel('$M_s$ (Nm)')
    plt.rcParams.update({'font.size': 14})


    plt.subplot(312)
    plt.loglog(M0[i],Dsub[i,0],'bo',markeredgecolor='black',markeredgewidth=0.2,alpha=0.5)
    plt.title('b)',loc='left')
    plt.ylabel('$T_S^0$ (s)')
    plt.grid(True,linewidth=0.25)
    plt.rcParams.update({'font.size': 14})


    plt.subplot(313)
    plt.semilogx(M0[i],Dsub[i,ibig]/Td[i]*100,'ro',markeredgecolor='black',markeredgewidth=0.2,alpha=0.3)
    plt.grid(True,linewidth=0.25)
    plt.ylim(0,100)
    plt.xlabel('$M_0$ (Nm)')
    plt.ylabel('Time of big subevent (s)')
    plt.title('c)',loc='left')
    plt.rcParams.update({'font.size': 14})
plt.savefig('FigureS4.pdf')


########## FIGURE S6 #############
fig=plt.figure(figsize=(8,10))
plt.subplot(211)
plt.hist(np.asarray(err[II]),100)
plt.title('a)',loc='left')
plt.rcParams.update({'font.size': 14})
plt.text(1.25,200,'median '+str(round(np.median(np.asarray(err[II]))*100)/100),fontsize=14)
plt.text(1.25,100,'std '+str(round(np.std(np.asarray(err[II]))*100)/100),fontsize=14)
plt.grid(True,linewidth=0.25)
plt.xlim(0.1,2)


data=np.load('allvar_triangle.npz') # load the triangle variables.
err=data['err']
Nsub=data['Nsub']
II=np.where(Nsub!=0)[0]
plt.subplot(212)
plt.hist(np.asarray(err[II]),100)
plt.title('b)',loc='left')
plt.rcParams.update({'font.size': 14})
plt.text(1.25,120,'median '+str(round(np.median(np.asarray(err[II]))*100)/100),fontsize=14)
plt.text(1.25,80,'std '+str(round(np.std(np.asarray(err[II]))*100)/100),fontsize=14)
plt.grid(True,linewidth=0.25)
plt.xlim(0.1,2)
plt.savefig('FigureS6.pdf')

plt.show()
