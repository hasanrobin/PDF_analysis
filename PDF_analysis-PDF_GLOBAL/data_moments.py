#!/usr/bin/python3

# Python package requirements
import os
import sys
import time
from datetime import datetime, date, timedelta
import numpy as np
import xarray as xr
import getopt
import json
import glob 
from scipy import stats
from math import isnan
from multiprocessing import Pool, cpu_count
from glob import glob
import pandas as pd 

# local requirements
from write_fitdistr import *
#from util_func import *
 
#---> Reading inputs for target year, month, variable 
year_start= sys.argv[1]
year_end =sys.argv[2]
start_mon= sys.argv[3]
end_mon=sys.argv[4]
input_var= sys.argv[5] 

yr_lst = []

for r in range(int(year_start),int(year_end)+1):
    for s in range (int(start_mon),int(end_mon)+1):
        yr_lst.append('%02i' %r + '%02i' %s)
print("yr-month=",yr_lst)

#--->reading 30yrs ECMWF files ----------

files_dir= "/data/oda/mg13420/model_wind/wind_91_20/"    # change path according to local directory 
   
file_list= sorted(list())
for (path, dirs, files) in os.walk(files_dir):
    file_list+=[file for file in files]
    #filelist= sorted(filelist)

target_files= list()
for f in file_list:
    for g in range(len(yr_lst)):
        if f.startswith((yr_lst[g])):
           target_files.append(f)
input_files= sorted(target_files)

i=0
files =[]
for i in range(len(input_files)):
   fl= files_dir+ input_files[i]
   files.append(fl)

print(len(files))
print(files[0])
print(files[-1])
ds0 =xr.open_dataset(files[-1])
latitude=  ds0.variables["lat"]
longitude= ds0.variables["lon"]
nlat=len(latitude)
nlon=len(longitude)
nday= len(files)
lsm=ds0.LSM.values
lsm[lsm[:,:,:]==1]=-500
lsm[lsm[:,:,:]==0]=1
lsm[:,0:31,0:56]=np.nan
lsm[:,0:20,90:146]=np.nan
lsm[:,0:54,265:357]=np.nan
lsm[0,60:80,300:350]=np.nan
lsm[0,100:110,100:120]=np.nan
lsm[0,118:128,340:350]=np.nan
lsm[0,130:137,310:340]=np.nan


d2m=xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32), dims=['time','lat','lon'])
t2m=xr.DataArray(np.zeros(((nday,nlat,nlon)),dtype=np.float32),  dims=['time','lat','lon'])
msl=xr.DataArray(np.zeros(((nday,nlat,nlon)),dtype=np.float32),  dims=['time','lat','lon'])
windu= xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32),dims=['time','lat','lon'])
windv= xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32),dims=['time','lat','lon'])
land = xr.DataArray(np.zeros((( nday, nlat, nlon)), dtype=np.float32), dims=['time','lat','lon'])
var_xr = xr.DataArray(np.zeros((( nday, nlat, nlon)), dtype=np.float32), dims=['time','lat','lon'])
for iday in  range(nday):
    file= xr.open_dataset(fls[iday])
    #d2m[iday,:,:]= np.mean(file.variables["D2M"][:,:,:],  axis=0)
    var_xr[iday,:,:]=  np.mean(file.variables[input_var][:,:,:], axis=0) # multiply with *0.01 for MSL-P variable 
    #msl[iday,:,:]=   np.mean(file.variables["MSL"][:,:,:], axis=0)*0.01
    #tcc[iday,:,:] =  np.mean(file.variables["TCC"][:,:,:], axis=0)
    #windu[iday,:,:]= np.mean(file.variables["U10M"][:,:,:], axis=0)
    #windv[iday,:,:]= np.mean(file.variables["V10M"][:,:,:], axis=0)
    #land[iday,:,:] =  np.mean( file.variables["LSM"][:,:,:],axis=0)
    iday+=1
#windm = np.sqrt((windu[:,:,:]**2 + windv[:,:,:]**2))


dtrend=False 
def anom_tseries(xarr,start_yr):
    dt=pd.date_range(start=start_yr, end='2020-12-31',freq='D')
    ds_xarr= xr.DataArray(xarr,dims=['time','lat','lon'], coords=dict(time=dt))
    ds_cl =ds_xarr.groupby('time.month').mean('time')
    ds_an= ds_xarr.groupby('time.month')-ds_cl
    return ds_an

ds_anom1 = anom_tseries(windu, "2006-01-01")
ds_anom2 = anom_tseries(windv, "2006-01-01")
#ds_anom  = np.sqrt((ds_anom1[:,:,:]**2 + ds_anom2[:,:,:]**2))

ds_anom= anom_tseries(var_xr,'2006-01-01')

atm_var = ds_anom
atm_len= len(atm_var[:,:,:])
crd_time=range(atm_len)
print("atm_var_shap",np.shape(atm_len)) 

# ----------------------------------------------------------------------
#       ***  main  funciton of the program.***
# ----------------------------------------------------------------------
if __name__=='__main__':

    time_start = time.time()
    date_start = datetime.now()

    print(" ------------------------------")
    print(" | MOMENTS OF PROB. DISTRIBUTION")
    print(" | python version="+sys.version)
    print(" | Start Date:", date_start)
    print(" ------------------------------")
    print("")
    log=None

    # Read script args
    idexp = ncores  = path_exp = None

    try:
        options, rem = getopt.getopt(sys.argv[1:], 'd:n:c:e:h', ['idexp=','ncores=','path_exp=','help'])

        print(options,rem)
        for opt, arg in options:
            if opt in ('-d', '--idexp'):
                idexp = arg
            if opt in ('-n', '--ncores'):
                ncores = int(arg)
       
            elif opt in ('-e', '--path_exp'):
                path_exp = arg
            elif opt in ('-h', '--help'):
                #printHelp(log)
                sys.exit(0)

    except getopt.GetoptError:
        print("ERROR! Wrong arguments!")
        printHelp(log)
        sys.exit(1)


    print(f" Script Args: ID-EXP = {idexp}")
    print(f'              NCORES = {ncores}')
    #print(f'              PATH-CONFIG = {path_config}')
    print(f'              PATH-EXP = {path_exp}')
    print("")

    #if (not idexp) or (not path_config) or (not path_exp):
    #    print("Wrong number of arguments!")
    #    printHelp(log)
    #    sys.exit(1)

    # If the experiment folder does not exist, create it
    #if not os.path.isdir(path_exp):
    #    print("ERROR! PATH "+ path_exp+" does not exisit!")
    #     sys.exit(1)

    #path_id_exp = path_exp + "/" + idexp
    # print("path_id_exp ="+path_id_exp)

    # If the experiment does not exist, exit
    #if not os.path.isdir(path_id_exp):
    #    print(f" ID-Experiment '{idexp}' DOES NOT EXISIT, make it!")
    #    sys.exit(1)

    # ------------------------------------------------------------------------------------
    print(f" Configuration Paramiters:")
    print(f" Time Series:   year_start = {year_start}")
    print(f'                year_end = {year_end}')
    #print(f'                month_start = {start_month}')
    #print(f'                month_end = {end_month}')
    #print(f'                day_start = {day_start}')
    #print(f'                day_end = {day_end}')
    #print(f" Datasets:      path_indata_base = {path_indata_base}")
    #print(f' Variables:     dim_lon = {dim_lon}')
    #print(f'                dim_lat = {dim_lat}')
    #print(f'                crd_lon = {crd_lon}')
    #print(f'                crd_lat = {crd_lat}')
    #print(f'                lvar_vector = {lvar_vector}')
    
    
    #if(lvar_vector):
    #    path_fitdata = path_id_exp+"/data/"+var_vectorMag
    #    file_momdata = "momData_"+var_vectorMag+".nc"
    #else:
    #    path_fitdata = path_id_exp+"/data/"+var_vector1
    #    file_momdata = "momData_"+var_vector1+".nc"

    #if not os.path.isdir(path_fitdata):
    #    print("ERROR! PATH "+ path_fitdata+" does not exisit!")
    #    sys.exit(1)
    
    var_name= str.lower(input_var)+"_anom" 
    file_dir= "/work/oda/mg13420/fitparm_exp_all/anom_06_20/"   # Change to local path/directory
    file_path= file_dir+var_name
   
    file_momdata= "momData_{}_06_20.nc".format(var_name) 
     
    #longitude = ds0.variables ["lon"]
    #latitude = ds0.variables["lat"][::-1]
    crd_lat= "latitude"
    crd_lon= "longitude"
 
    nlat=len(latitude)
    nlon= len(longitude)
    var_vector1 = "wndm"
    var2 =  "land"
    nday= len(crd_time)
    #dim_lat= "latitude" 
    #dim_lon= "longitude"

    #-- read variable latitude and longitude arrays
    #with xr.open_dataset(path, engine= "netcdf4") as ds0:
        # print(ds0)
    #    longitude = ds0[crd_lon]
    #    latitude  = ds0[crd_lat]
    #    times     = ds0[crd_time]
    #    nlon  = len(longitude)
    #    nlat  = len(latitude)
    #    ntime = len(times)

    #    if (len(ds0.dims) == 3):
    #        ndepth = 1
    #    elif (len(ds0.dims) == 4):
    #        depth  = ds0[crd_depth]
    #        ndepth = len(depth)
    ndepth=1
   
    print(f" Data grid size:")
    print(f"    nlon   = {nlon}")
    print(f"    nlat   = {nlat}")
    print(f"    ntime_day  = {nday}")
   
    print("")
    
    field = xr.DataArray(np.zeros((nday,nlat,nlon)),dims=['time','lat','lon'])
    #lsm = np.zeros((nday, nlat,nlon)) # ,dims=['time','depth','lat','lon'])
    
    field[:,:,:] = pdf_var[:,:,:] # * lsm[0,:,:]     #---> multiply with lsm for Land-Sea mask  
    
   #Now assign the results to your array
    if (len(ds0.dims) == 3):
        mean     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        variance      = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        skew     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        kurtosis = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
    elif (len(ds0.dims) == 4):
        mean     = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
        variance      = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
        skew     = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
        kurtosis = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
    else:
        print(f" Error! Dimension not correct")

    #lsmask =np.ma.masked_where(lsm==0, lsm, copy=False )
    mean.attrs['long_name']     = 'mean'
    variance.attrs['long_name'] = 'variance'
    skew.attrs['long_name']     = 'skew'
    kurtosis.attrs['long_name'] = 'kurtosis'
     
    if (len(ds0.dims) == 3):
       for ij in range(1):
           for ik in range(nlat):
               for il in range(nlon):
                    if (lsm[ij,ik,il]==-500):continue 
                    #if not isnan(field[ij,ik,il]):
 
                    mean[ik,il]     = np.mean(field[:,ik,il], axis = 0)
                    variance[ik,il]   = np.var(field[:,ik,il], axis = 0)
                    skew[ik,il]     = stats.skew(field[:,ik,il], axis=0)
                    kurtosis[ik,il] = stats.kurtosis(field[:,ik,il],fisher=True,bias=True, axis=0) 

                    write_moments_2D(file_path,file_momdata,longitude,latitude,mean,variance,skew, kurtosis)


    time_end = time.time()
    elapsed_time = (time_end - time_start)
    print(file_path+'/'+file_momdata) 
    print(" elapsed time {}".format(elapsed_time)+"sec = {}".format(elapsed_time/60.)+"min")
