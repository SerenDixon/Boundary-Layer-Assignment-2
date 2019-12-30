# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:38:22 2019

@author: 25006773
"""

from netCDF4 import Dataset as ncfile
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Give the location of the MJO file and open it
infile='/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/MJO_rmm1_rmm2.jan-dec_dmeans_ts.1979-2018 (1) - Copy.nc'
nc_infile = ncfile(infile)

# Read in variables
mjo_phase = nc_infile['phase_ts'][:]
mjo_amp= nc_infile['amplitude_ts'][:]
rmm1= nc_infile['rmm1_ts'][:]
rmm2= nc_infile['rmm2_ts'][:]
#Read in time variable to get list of dates
time = nc_infile['time'][:]

#create empty MJO date lists
phase0_dates=[]
phase1_dates=[]
phase2_dates=[]
phase3_dates=[]
phase4_dates=[]
phase5_dates=[]
phase6_dates=[]
phase7_dates=[]
phase8_dates=[]

#to sort out dates of each phase into lists, selecting dates only with an amplitude larger than one
for i in time:
    if mjo_phase[i]==0 and mjo_amp[i]>=1.0:
        phase0_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==1 and mjo_amp[i]>=1.0:
        phase1_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==2 and mjo_amp[i]>=1.0:
        phase2_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==3 and mjo_amp[i]>=1.0:
        phase3_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==4 and mjo_amp[i]>=1.0:
        phase4_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==5 and mjo_amp[i]>=1.0:
        phase5_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==6 and mjo_amp[i]>=1.0:
        phase6_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==7 and mjo_amp[i]>=1.0:
        phase7_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_phase[i]==8 and mjo_amp[i]>=1.0:
        phase8_dates.append(dt.timedelta(days=float(i))+dt.datetime(year=1979,month=1,day=1))
        i==i+1
    elif mjo_amp[i]<1.0:
        i=i+1

#find dates for phases 4 and 8        
print (phase4_dates)
print (phase8_dates)

#dates selected from those lists and the relevant data was downloaded from AMDAR

#read in temp and pressure data for phases 4 and 8
phase4_1=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p4_0404_csv.csv')
phase4_2=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p4_0832_csv.csv')
phase4_3=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p4_midday_csv.csv')
phase4_4=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p4_1631_csv.csv')
phase4_5=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p4_2346_csv.csv')

phase8_1=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p8_0419_csv.csv')
phase8_2=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p8_0840_csv.csv')
phase8_3=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p8_midday_csv.csv')
phase8_4=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p8_1649_csv.csv')
phase8_5=pd.read_csv('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Data_Sing_p8_2253_csv.csv')

#split out data from columns

#dry temp

#phase 4
t4_1=phase4_1.loc[:,'t']
t4_2=phase4_2.loc[:,'t']
t4_3=phase4_3.loc[:,'t']
t4_4=phase4_4.loc[:,'t']
t4_5=phase4_5.loc[:,'t']

#phase 8
t8_1=phase8_1.loc[:,'t']
t8_2=phase8_2.loc[:,'t']
t8_3=phase8_3.loc[:,'t']
t8_4=phase8_4.loc[:,'t']
t8_5=phase8_5.loc[:,'t']

#dewpoint temp

#phase 4
td4_1=phase4_1.loc[:,'td']
td4_2=phase4_2.loc[:,'td']
td4_3=phase4_3.loc[:,'td']
td4_4=phase4_4.loc[:,'td']
td4_5=phase4_5.loc[:,'td']

#phase 8
td8_1=phase8_1.loc[:,'td']
td8_2=phase8_2.loc[:,'td']
td8_3=phase8_3.loc[:,'td']
td8_4=phase8_4.loc[:,'td']
td8_5=phase8_5.loc[:,'td']


#height in ft

#phase 4
h4_1=phase4_1.loc[:,'Plane altitude (ft)']
h4_2=phase4_2.loc[:,'Plane altitude (ft)']
h4_3=phase4_3.loc[:,'Plane altitude (ft)']
h4_4=phase4_4.loc[:,'Plane altitude (ft)']
h4_5=phase4_5.loc[:,'Plane altitude (ft)']

#phase 8
h8_1=phase8_1.loc[:,'Plane altitude (ft)']
h8_2=phase8_2.loc[:,'Plane altitude (ft)']
h8_3=phase8_3.loc[:,'Plane altitude (ft)']
h8_4=phase8_4.loc[:,'Plane altitude (ft)']
h8_5=phase8_5.loc[:,'Plane altitude (ft)']

#pressure

#phase 4
p4_1=phase4_1.loc[:,'mb']
p4_2=phase4_2.loc[:,'mb']
p4_3=phase4_3.loc[:,'mb']
p4_4=phase4_4.loc[:,'mb']
p4_5=phase4_5.loc[:,'mb']

#phase 8
p8_1=phase8_1.loc[:,'mb']
p8_2=phase8_2.loc[:,'mb']
p8_3=phase8_3.loc[:,'mb']
p8_4=phase8_4.loc[:,'mb']
p8_5=phase8_5.loc[:,'mb']

#CALCULATE POTENTIAL TEMP

#first need surface pressure from each time for each date (in mb)

#phase 4 (18/1/18)
sp_4_1=1007.4#0404
sp_4_2=1008.5#0832
sp_4_3=1008.5#1239
sp_4_4=1006.4#1631
sp_4_5=1009.5#2346

#phase 8 (27/12/17)
sp_8_1=1008.5#0419
sp_8_2=1009.5#0840
sp_8_3=1010.5#1118
sp_8_4=1006.4#1649
sp_8_5=1009.5#2253

#also need R - the gas constant for dry air
#and cp - the gas constant for dry air at constant volume
R=287 #J K-1 kg-1
cp=1004 #J/kg.K

#now need to convert temps into kelvin

#dry temps

#phase4 (times 1 to 5)

t4_1_k=[] 
i=0
for i in t4_1:
    a=i+273.15
    t4_1_k.append(a)
    i+1
    
t4_2_k=[] 
i=0
for i in t4_2:
    a=i+273.15
    t4_2_k.append(a)
    i+1
    
t4_3_k=[] 
i=0
for i in t4_3:
    a=i+273.15
    t4_3_k.append(a)
    i+1
    
t4_4_k=[] 
i=0
for i in t4_4:
    a=i+273.15
    t4_4_k.append(a)
    i+1
    
t4_5_k=[] 
i=0
for i in t4_5:
    a=i+273.15
    t4_5_k.append(a)
    i+1
    
#phase8 (times 1 to 5)

t8_1_k=[] 
i=0
for i in t8_1:
    a=i+273.15
    t8_1_k.append(a)
    i+1
    
t8_2_k=[] 
i=0
for i in t8_2:
    a=i+273.15
    t8_2_k.append(a)
    i+1
    
t8_3_k=[] 
i=0
for i in t8_3:
    a=i+273.15
    t8_3_k.append(a)
    i+1
    
t8_4_k=[] 
i=0
for i in t8_4:
    a=i+273.15
    t8_4_k.append(a)
    i+1
    
t8_5_k=[] 
i=0
for i in t8_5:
    a=i+273.15
    t8_5_k.append(a)
    i+1
    
#dewpoint temps
    
#phase4 (times 1 to 5)
    
td4_1_k=[] 
i=0
for i in td4_1:
    a=i+273.15
    td4_1_k.append(a)
    i+1
    
td4_2_k=[] 
i=0
for i in td4_2:
    a=i+273.15
    td4_2_k.append(a)
    i+1
    
td4_3_k=[] 
i=0
for i in td4_3:
    a=i+273.15
    td4_3_k.append(a)
    i+1
    
td4_4_k=[] 
i=0
for i in td4_4:
    a=i+273.15
    td4_4_k.append(a)
    i+1
    
td4_5_k=[] 
i=0
for i in td4_5:
    a=i+273.15
    td4_5_k.append(a)
    i+1
    
#phase8 (times 1 to 5)
    
td8_1_k=[] 
i=0
for i in td8_1:
    a=i+273.15
    td8_1_k.append(a)
    i+1
    
td8_2_k=[] 
i=0
for i in td8_2:
    a=i+273.15
    td8_2_k.append(a)
    i+1
    
td8_3_k=[] 
i=0
for i in td8_3:
    a=i+273.15
    td8_3_k.append(a)
    i+1
    
td8_4_k=[] 
i=0
for i in td8_4:
    a=i+273.15
    td8_4_k.append(a)
    i+1
    
td8_5_k=[] 
i=0
for i in td8_5:
    a=i+273.15
    td8_5_k.append(a)
    i+1

#create lists of pressures 

#phase4 (times 1 to 5)
    
P4_1=[]
i=0
for i in p4_1:
    P4_1.append(i)
    i+1
    
P4_2=[]
i=0
for i in p4_2:
    P4_2.append(i)
    i+1
    
P4_3=[]
i=0
for i in p4_3:
    P4_3.append(i)
    i+1
    
P4_4=[]
i=0
for i in p4_4:
    P4_4.append(i)
    i+1
    
P4_5=[]
i=0
for i in p4_5:
    P4_5.append(i)
    i+1
    
#phase8 (times 1 to 5)
    
P8_1=[]
i=0
for i in p8_1:
    P8_1.append(i)
    i+1
    
P8_2=[]
i=0
for i in p8_2:
    P8_2.append(i)
    i+1
    
P8_3=[]
i=0
for i in p8_3:
    P8_3.append(i)
    i+1
    
P8_4=[]
i=0
for i in p8_4:
    P8_4.append(i)
    i+1
    
P8_5=[]
i=0
for i in p8_5:
    P8_5.append(i)
    i+1
    
#now calculate pot temp
 
#phase4 (times 1 to 5)
    
pt4_1=[]
i=0.0
b=0
for i in t4_1_k and P4_1:
    a=t4_1_k[b]*((sp_4_1/P4_1[b])**(R/cp))
    pt4_1.append(a)
    b=b+1
    i=i+1
    
pt4_2=[]
i=0.0
b=0
for i in t4_2_k and P4_2:
    a=t4_2_k[b]*((sp_4_2/P4_2[b])**(R/cp))
    pt4_2.append(a)
    b=b+1
    i=i+1
    
pt4_3=[]
i=0.0
b=0
for i in t4_3_k and P4_3:
    a=t4_3_k[b]*((sp_4_3/P4_3[b])**(R/cp))
    pt4_3.append(a)
    b=b+1
    i=i+1
    
pt4_4=[]
i=0.0
b=0
for i in t4_4_k and P4_4:
    a=t4_4_k[b]*((sp_4_4/P4_4[b])**(R/cp))
    pt4_4.append(a)
    b=b+1
    i=i+1
    
pt4_5=[]
i=0.0
b=0
for i in t4_5_k and P4_5:
    a=t4_5_k[b]*((sp_4_5/P4_5[b])**(R/cp))
    pt4_5.append(a)
    b=b+1
    i=i+1

#phase8 (times 1 to 5)
    
pt8_1=[]
i=0.0
b=0
for i in t8_1_k and P8_1:
    a=t8_1_k[b]*((sp_8_1/P8_1[b])**(R/cp))
    pt8_1.append(a)
    b=b+1
    i=i+1
    
pt8_2=[]
i=0.0
b=0
for i in t8_2_k and P8_2:
    a=t8_2_k[b]*((sp_8_2/P8_2[b])**(R/cp))
    pt8_2.append(a)
    b=b+1
    i=i+1
    
pt8_3=[]
i=0.0
b=0
for i in t8_3_k and P8_3:
    a=t8_3_k[b]*((sp_8_3/P8_3[b])**(R/cp))
    pt8_3.append(a)
    b=b+1
    i=i+1
    
pt8_4=[]
i=0.0
b=0
for i in t8_4_k and P8_4:
    a=t8_4_k[b]*((sp_8_4/P8_4[b])**(R/cp))
    pt8_4.append(a)
    b=b+1
    i=i+1
    
pt8_5=[]
i=0.0
b=0
for i in t8_5_k and P8_5:
    a=t8_5_k[b]*((sp_8_5/P8_5[b])**(R/cp))
    pt8_5.append(a)
    b=b+1
    i=i+1

#convert heights from feet to meters

#phase4 (times 1 to 5)
    
h4_1_m=[] 
i=0 
for i in h4_1:
    a=i*0.3048
    h4_1_m.append(a)
    i+1
    
h4_2_m=[] 
i=0 
for i in h4_2:
    a=i*0.3048
    h4_2_m.append(a)
    i+1

h4_3_m=[] 
i=0 
for i in h4_3:
    a=i*0.3048
    h4_3_m.append(a)
    i+1
    
h4_4_m=[] 
i=0 
for i in h4_4:
    a=i*0.3048
    h4_4_m.append(a)
    i+1
    
h4_5_m=[] 
i=0 
for i in h4_5:
    a=i*0.3048
    h4_5_m.append(a)
    i+1
    
#phase8 (times 1 to 5)
    
h8_1_m=[] 
i=0 
for i in h8_1:
    a=i*0.3048
    h8_1_m.append(a)
    i+1
    
h8_2_m=[] 
i=0 
for i in h8_2:
    a=i*0.3048
    h8_2_m.append(a)
    i+1

h8_3_m=[] 
i=0 
for i in h8_3:
    a=i*0.3048
    h8_3_m.append(a)
    i+1
    
h8_4_m=[] 
i=0 
for i in h8_4:
    a=i*0.3048
    h8_4_m.append(a)
    i+1
    
h8_5_m=[] 
i=0 
for i in h8_5:
    a=i*0.3048
    h8_5_m.append(a)
    i+1
    
#plot pot temp against height 

#phase 4 time 1
plt.figure()
plt.plot(pt4_1,h4_1_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_1')
plt.show() 

#phase 4 time 2
plt.figure()
plt.plot(pt4_2,h4_2_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_2')
plt.show() 

#phase 4 time 3
plt.figure()
plt.plot(pt4_3,h4_3_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_3')
plt.show() 

#phase 4 time 4
plt.figure()
plt.plot(pt4_4,h4_4_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_4')
plt.show() 

#phase 4 time 5
plt.figure()
plt.plot(pt4_5,h4_5_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_5')
plt.show() 

#phase 4 all times
plt.figure()
plt.plot(pt4_1,h4_1_m, marker='o')
plt.plot(pt4_2,h4_2_m, marker='o')
plt.plot(pt4_3,h4_3_m, marker='o')
plt.plot(pt4_4,h4_4_m, marker='o')
plt.plot(pt4_5,h4_5_m, marker='o')
plt.legend('12345')
plt.ylim(0,3000)
plt.ylabel('Height (m)')
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_4_all')
plt.show()

#phase 8 time 1
plt.figure()
plt.plot(pt8_1,h8_1_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_1')
plt.show() 

#phase 8 time 2
plt.figure()
plt.plot(pt8_2,h8_2_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_2')
plt.show()

#phase 8 time 3
plt.figure()
plt.plot(pt8_3,h8_3_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_3')
plt.show()

#phase 8 time 4
plt.figure()
plt.plot(pt8_4,h8_4_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_4')
plt.show()

#phase 8 time 5
plt.figure()
plt.plot(pt8_5,h8_5_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_5')
plt.show()

#phase 8 all times
plt.figure()
plt.plot(pt8_1,h8_1_m, marker='o')
plt.plot(pt8_2,h8_2_m, marker='o')
plt.plot(pt8_3,h8_3_m, marker='o')
plt.plot(pt8_4,h8_4_m, marker='o')
plt.plot(pt8_5,h8_5_m, marker='o')
plt.legend('12345')
plt.ylim(0,3000)
plt.ylabel('Height (m)')
plt.xlabel('Potential temperature (K)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/PT_8_all')
plt.show()

#TO FIND RELATIVE HUMIDITY PROFILES

#first find saturation vapour pressure of air at t:

es0=6.11#mb (Saturation vapour pressure at 273K)
L=2500000#Jkg-1 (Latent heat of vaporisation)
rv=461#Jkg-1K-1 (Gas constant for water vapour)
t0=273#K

#phase4 (times 1 to 5)

es4_1_t=[] 
i=0
for i in t4_1_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_1_t.append(a)
    i+1
    
es4_2_t=[] 
i=0
for i in t4_2_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_2_t.append(a)
    i+1
    
es4_3_t=[] 
i=0
for i in t4_3_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_3_t.append(a)
    i+1
    
es4_4_t=[] 
i=0
for i in t4_4_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_4_t.append(a)
    i+1
    
es4_5_t=[] 
i=0
for i in t4_5_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_5_t.append(a)
    i+1
    
#phase8 (times 1 to 5)

es8_1_t=[] 
i=0
for i in t8_1_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_1_t.append(a)
    i+1
    
es8_2_t=[] 
i=0
for i in t8_2_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_2_t.append(a)
    i+1
    
es8_3_t=[] 
i=0
for i in t8_3_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_3_t.append(a)
    i+1
    
es8_4_t=[] 
i=0
for i in t8_4_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_4_t.append(a)
    i+1
    
es8_5_t=[] 
i=0
for i in t8_5_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_5_t.append(a)
    i+1
    
#now find saturation vapour pressure of air at td:

#phase4 (times 1 to 5)
    
es4_1_td=[] 
i=0
for i in td4_1_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_1_td.append(a)
    i+1
    
es4_2_td=[] 
i=0
for i in td4_2_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_2_td.append(a)
    i+1
    
es4_3_td=[] 
i=0
for i in td4_3_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_3_td.append(a)
    i+1
    
es4_4_td=[] 
i=0
for i in td4_4_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_4_td.append(a)
    i+1
    
es4_5_td=[] 
i=0
for i in td4_5_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es4_5_td.append(a)
    i+1
    
#phase8 (times 1 to 5)
    
es8_1_td=[] 
i=0
for i in td8_1_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_1_td.append(a)
    i+1
    
es8_2_td=[] 
i=0
for i in td8_2_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_2_td.append(a)
    i+1
    
es8_3_td=[] 
i=0
for i in td8_3_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_3_td.append(a)
    i+1
    
es8_4_td=[] 
i=0
for i in td8_4_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_4_td.append(a)
    i+1
    
es8_5_td=[] 
i=0
for i in td8_5_k:
    a=es0*(np.exp((L/rv)*((1/273)-(1/i))))
    es8_5_td.append(a)
    i+1
    
#now find relative humidity

#phase4 (times 1 to 5)
    
rh4_1=[] 
i=0
b=0
for i in es4_1_td and es4_1_t:
    a=es4_1_td[b]/es4_1_t[b]
    a=a*100
    rh4_1.append(a)
    b=b+1
    i+1
    
rh4_2=[] 
i=0
b=0
for i in es4_2_td and es4_2_t:
    a=es4_2_td[b]/es4_2_t[b]
    a=a*100
    rh4_2.append(a)
    b=b+1
    i+1
    
rh4_3=[] 
i=0
b=0
for i in es4_3_td and es4_3_t:
    a=es4_3_td[b]/es4_3_t[b]
    a=a*100
    rh4_3.append(a)
    b=b+1
    i+1
    
rh4_4=[] 
i=0
b=0
for i in es4_4_td and es4_4_t:
    a=es4_4_td[b]/es4_4_t[b]
    a=a*100
    rh4_4.append(a)
    b=b+1
    i+1
    
rh4_5=[] 
i=0
b=0
for i in es4_5_td and es4_5_t:
    a=es4_5_td[b]/es4_5_t[b]
    a=a*100
    rh4_5.append(a)
    b=b+1
    i+1
    
#phase8 (times 1 to 5)
    
rh8_1=[] 
i=0
b=0
for i in es8_1_td and es8_1_t:
    a=es8_1_td[b]/es8_1_t[b]
    a=a*100
    rh8_1.append(a)
    b=b+1
    i+1
    
rh8_2=[] 
i=0
b=0
for i in es8_2_td and es8_2_t:
    a=es8_2_td[b]/es8_2_t[b]
    a=a*100
    rh8_2.append(a)
    b=b+1
    i+1
    
rh8_3=[] 
i=0
b=0
for i in es8_3_td and es8_3_t:
    a=es8_3_td[b]/es8_3_t[b]
    a=a*100
    rh8_3.append(a)
    b=b+1
    i+1
    
rh8_4=[] 
i=0
b=0
for i in es8_4_td and es8_4_t:
    a=es8_4_td[b]/es8_4_t[b]
    a=a*100
    rh8_4.append(a)
    b=b+1
    i+1
    
rh8_5=[] 
i=0
b=0
for i in es8_5_td and es8_5_t:
    a=es8_5_td[b]/es8_5_t[b]
    a=a*100
    rh8_5.append(a)
    b=b+1
    i+1
    
#plot relative humudity against height 

#phase 4 time 1
plt.figure()
plt.plot(rh4_1,h4_1_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_1')
plt.show() 

#phase 4 time 2
plt.figure()
plt.plot(rh4_2,h4_2_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_2')
plt.show()   

#phase 4 time 3
plt.figure()
plt.plot(rh4_3,h4_3_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_3')
plt.show()  

#phase 4 time 4
plt.figure()
plt.plot(rh4_4,h4_4_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_4')
plt.show()  

#phase 4 time 5
plt.figure()
plt.plot(rh4_5,h4_5_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_5')
plt.show()  

#phase 4 all times
plt.figure()
plt.plot(rh4_1,h4_1_m, marker='o')
plt.plot(rh4_2,h4_2_m, marker='o')
plt.plot(rh4_3,h4_3_m, marker='o')
plt.plot(rh4_4,h4_4_m, marker='o')
plt.plot(rh4_5,h4_5_m, marker='o')
plt.legend('12345')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative Humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_4_all')
plt.show()

#phase 8 time 1
plt.figure()
plt.plot(rh8_1,h8_1_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_1')
plt.show() 

#phase 8 time 2
plt.figure()
plt.plot(rh8_2,h8_2_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_2')
plt.show() 

#phase 8 time 3
plt.figure()
plt.plot(rh8_3,h8_3_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_3')
plt.show()

#phase 8 time 4
plt.figure()
plt.plot(rh8_4,h8_4_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_4')
plt.show() 

#phase 8 time 5
plt.figure()
plt.plot(rh8_5,h8_5_m, marker='o')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_5')
plt.show() 

#phase 8 all times
plt.figure()
plt.plot(rh8_1,h8_1_m, marker='o')
plt.plot(rh8_2,h8_2_m, marker='o')
plt.plot(rh8_3,h8_3_m, marker='o')
plt.plot(rh8_4,h8_4_m, marker='o')
plt.plot(rh8_5,h8_5_m, marker='o')
plt.legend('12345')
plt.ylabel('Height (m)')
plt.ylim(0,3000)
plt.xlabel('Relative Humidity (%)')
plt.savefig('/Users/seren/Documents/University/Year 3/MT37J (Boundary Layer Met)/Assignment 2/Graphs/RH_8_all')
plt.show()
