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
in_yr= sys.argv[1]
end_yr =sys.argv[2]
in_mon= sys.argv[3]
end_mon=sys.argv[4]
input_var= sys.argv[5] 

start_day=date(int(in_yr),int(in_mon),1)
start_date= start_day.strftime("%Y-%m-%d")
end_day=date(int(end_yr),int(end_mon),31)
print(start_date)
end_date= end_day.strftime("%Y-%m-%d")
print(end_date)

yr_lst = []

for r in range(int(in_yr),int(end_yr)+1):
    for s in range (int(in_mon),int(end_mon)+1):
        yr_lst.append('%02i' %r + '%02i' %s)
print("yr-month=",yr_lst)

#--->reading 30yrs ECMWF files ----------

files_dir= "/data/cmcc/mg13420/sst_data/sst_IBI_0.05deg/"    # change path according to local directory 
#files_dir="/data/cmcc/mg13420/sst_data/obs_sst_IBI_0.05deg/"

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
nday= len(files)
print("ds_shape",len(ds0.dims)) 
sst_comb=xr.open_mfdataset(files,combine='nested', parallel=False , concat_dim= ['time'], )
#sst= sst_comb['analysed_sst'][:,:,:]  
sst= sst_comb['thetao'][:,0,:,:]
#sst = xr.DataArray(sst, dims=['time','latitude','longitude'], coords=[atm_comb['time'][:8], latitude, longitude]) 
sst = xr.DataArray(sst, dims=['time','lat','lon'], coords=[sst_comb['time'], latitude, longitude]) 
#sst =np.nan_to_num(sst, nan=-500)
dtrend=False 
#-- processing for anomaly time series          
def anom_tseries(xarr,start_date, end_date):
    dt= pd.date_range(start=start_date, end=end_date, freq='D') # freq=daily 
    ds_xarr= xr.DataArray(xarr,dims=['time','lat','lon'], coords=dict(time=dt))
    ds_cl =ds_xarr.groupby('time.dayofyear').mean('time') 
    ds_an= ds_xarr.groupby('time.dayofyear')-ds_cl
    return ds_an

anom_var = anom_tseries(sst, start_date, end_date)
print("Anomaly shape=", anom_var[0].values)

sst_var = sst
sst_len= len(sst_var[:,:,:])
crd_time=range(sst_len)
print("atm_var_shap",np.shape(sst_len)) 

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
    # idexp = ncores  = path_exp = None

    # try:
    #     options, rem = getopt.getopt(sys.argv[1:], 'd:n:c:e:h', ['idexp=','ncores=','path_exp=','help'])

    #     print(options,rem)
    #     for opt, arg in options:
    #         if opt in ('-d', '--idexp'):
    #             idexp = arg
    #         if opt in ('-n', '--ncores'):
    #             ncores = int(arg)
       
    #         elif opt in ('-e', '--path_exp'):
    #             path_exp = arg
    #         elif opt in ('-h', '--help'):
    #             #printHelp(log)
    #             sys.exit(0)

    # except getopt.GetoptError:
    #     print("ERROR! Wrong arguments!")
    #     printHelp(log)
    #     sys.exit(1)


    # print(f" Script Args: ID-EXP = {idexp}")
    # print(f'              NCORES = {ncores}')
    # #print(f'              PATH-CONFIG = {path_config}')
    # print(f'              PATH-EXP = {path_exp}')
    # print("")

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

   

    # ------------------------------------------------------------------------------------
    print(f" Configuration Paramiters:")
    print(f" Time Series:   year_start = {in_yr}")
    print(f'                year_end = {end_yr}')
    print(f'                month_start = {in_mon}')
    print(f'                month_end = {end_mon}')
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
    
    var_name= str.lower(input_var)# +"_anom" 
    file_dir= "/work/cmcc/mg13420/plot_exercises/sst_pdf-exps/"   # Change to local path/directory
    file_path= file_dir+var_name
   
    #file_momdata= "momData_{}_2018_2023.nc".format(var_name) 
    file_momdata= "momData_{}_2022_2023.nc".format(var_name) 

    
    crd_lat= "latitude"
    crd_lon= "longitude"
    nlat=len(latitude)
    nlon= len(longitude)
    nday= len(crd_time)
    

    # if (len(ds0.dims) == 3):
    #     ndepth = 1
    # elif (len(ds0.dims) == 4):
    #     depth  = ds0[crd_depth]
    # #        ndepth = len(depth)
    ndepth=1
   
    print(f" Data grid size:")
    print(f"    nlon   = {nlon}")
    print(f"    nlat   = {nlat}")
    print(f"    ntime_day  = {nday}")
   
    print("")
    
    field = xr.DataArray(np.zeros((nday,nlat,nlon)),dims=['time','lat','lon'])
    
    field[:,:,:] = anom_var[:,:,:] # * lsm[0,:,:]     #---> multiply with lsm for Land-Sea mask  
    
   #Now assign the results to your array
    if (len(ds0.dims) == 4):
        mean     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        variance      = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        skew     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        kurtosis = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
    # elif (len(ds0.dims) == 4):
    #     mean     = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
    #     variance      = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
    #     skew     = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
    #     kurtosis = xr.DataArray(np.zeros((ndepth,nlat,nlon)), coords=[depth, latitude, longitude], dims=['depth','latitude','longitude'])
    # else:
    #     print(f" Error! Dimension not correct")

    #lsmask =np.ma.masked_where(lsm==0, lsm, copy=False )
    mean.attrs['long_name']     = 'mean'
    variance.attrs['long_name'] = 'variance'
    skew.attrs['long_name']     = 'skew'
    kurtosis.attrs['long_name'] = 'kurtosis'
     
    if (len(ds0.dims) == 4):
       # for ij in range(1):
       for ik in range(nlat):
           for il in range(nlon):
                if np.isfinite(field[0,ik,il]):#==-500):
                    #continue 
                #if not isnan(field[ij,ik,il]):
                    mean[ik,il]     = np.mean(field[:,ik,il], axis = 0)
                    variance[ik,il]   = np.var(field[:,ik,il], axis = 0)
                    skew[ik,il]     = stats.skew(field[:,ik,il], axis=0)
                    kurtosis[ik,il] = stats.kurtosis(field[:,ik,il],fisher=True,bias=True, axis=0) 

    write_moments_2D(file_path,file_momdata,longitude,latitude, crd_lat, crd_lon,mean,variance,skew, kurtosis)


    time_end = time.time()
    elapsed_time = (time_end - time_start)
    print(file_path+'/'+file_momdata) 
    print(" elapsed time {}".format(elapsed_time)+"sec = {}".format(elapsed_time/60.)+" min")
