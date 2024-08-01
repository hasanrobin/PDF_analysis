#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from netCDF4 import Dataset
import numpy as np
from glob import glob
import xarray as xr 
import time
import logging 
import shutil 
from scipy.stats import exponweib
from datetime import datetime, date, timedelta
from multiprocessing import Pool, cpu_count
from scipy import stats
from math import isnan 
import  pandas as pd

# ----------------------------------------------------------------------
#                ***  Parallel processing_function  ***
# Purpose: divide the input time series based on spatial domain size
# ----------------------------------------------------------------------

def worker_func(mylist):

    process_id = os.getpid()
    start_time = datetime.now()
    log.info(" proc-ID "+str(process_id)+" started at {}".format(start_time))
    print(" proc-ID "+str(process_id)+" started at {}".format(start_time))

    nlon_sub = len(mylist[0])
    crd_lon     = mylist[1]
    crd_lat     = mylist[2]
    crd_time    = mylist[3] 
    #files =       mylist[4]
    #var_vector1 = mylist[5]

    ds0 = xr.open_dataset(files[0],decode_times=False )
    #print(crd_time)
    ndays    =   len(crd_time)
    mlat     =  len(crd_lat)
    mlon     =  len(crd_lon)
    
# -----------------------------------------------------
# Define/upload the variable 'field'/values to be fitted
# ----------------------------------------------------- 
    
    field = xr.DataArray(np.zeros(((ndays, mlat,nlon_sub)), dtype=np.float64),dims=['time', 'lat','lon']) 
    field[:,:,:] = pdf_var.isel(longitude=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[:,:,:] # Use this if we broadcast input values from the main function section
    #lsm =   xr.DataArray(np.zeros(((4,mlat,nlon_sub)), dtype=np.float32), dims=['time', 'lat','lon'])
    #lsm[:,:,:] = land.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[:,:,:]
                    
    # if (len(ds0.dims) == 3):
    #     if(lvar_vector):
    #         for iday in range(len(files)):
    #             ds = xr.open_dataset(files[iday])
    #             vector1 = np.mean(ds.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[var_vector1], axis=0)
    #             vector2 = np.mean(ds.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[var_vector2], axis=0 )
    #             vector1 = ds.variables[var_vector1].values
    #             vector2 = ds.variables[var_vector2].values
    #             field[iday,:,:] = np.sqrt(vector1**2+vector2**2)
    #             print(np.shape(field))
    #             iday+= 1
    #     else:
    #         for iday in range(ndays):
    #             ds = xr.open_dataset(files[0] )
    #             field[iday,:,:] = ds.isel(longitude=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[var_vector1]# multiply with *0.01 for MSLP
    #             print(f"Shape of field[{iday},:,:]: {np.shape(field[iday,:,:])}")
    #             iday+=1
            
    
# Define the FIT parameters, PDF  statistics
# -----------------------------------------
    
    pdf_name= getattr(stats,pdf_fit)
    # check the number of parameters of the PDF
    nparams = pdf_name.numargs
    print(f"Number of parameters for {pdf_fit} is {nparams}")
    

    pdf_params = xr.DataArray(np.zeros((nparams, mlat,nlon_sub)), dtype=np.float64)
    #pdf_params[:,:,:] = np.nan
    print(np.shape(pdf_params))
    #  fit the PDF to the data depending on the number of parameters 
    if (nparams == 2):
        for iy in range(mlat):
            for ix in range(nlon_sub):
                #if not isnan(field[0,iy,ix]):
                pdf_params[:,iy,ix] = list(pdf_name.fit(field[:,iy,ix]))
                
    if (nparams == 3):
        for iy in range(mlat):
            for ix in range(nlon_sub):
                #if not isnan(field[0,iy,ix]):
                pdf_params[:,iy,ix] = list(pdf_name.fit(field[:,iy,ix],floc=0))

    # if (pdf_name == "weib"):
    #     for iz in range(ndepth):
    #         for iy in range(mlat):
    #             for ix in range(nlon_sub):
    #                 #if (lsm[iz, iy, ix]==-500) :continue
    #                 #if not isnan(field[0,0,iy,ix]):
    #                 pdf_params[:,iy,ix] = list(stats.weibull_min.fit(field[:,iy,ix],loc=0.01))
    
    # #---This skewnormal PDF function is adapted with the Nelder-med algorithm, if we have convergent issue in shape parameter distribution  
    # if (pdf_name == "skew"):
    #         for iy in range(mlat):
    #             for ix in range(nlon_sub):
    #             # print(f'ix={ix}')
    #                 #if (lsm[0,iy,ix]==-500) :continue
    #                 #if not isnan(field[iz,iy,ix]):
    #                 n=2    # We can try with different number set 
    #                 par_min=  [-6,field.min(),5] 	
    #                 par_max=  [1, field.max(),10] 
    #                 rand_parm =np.random.uniform(low=par_min, high =par_max, size =(n,3))
    #                 #rand_parm= np.round(rand_parm,2
    #                 field_min=10**10
    #                 for ipar in range(len(rand_parm)):
    #                     results = list(stats.skewnorm.fit(field[:,iy,ix],rand_parm[ipar][0],loc=rand_parm[ipar][1],scale=rand_parm[ipar][2] ))
    #                     parfit = results[0]
    #                     if (results[1]< field_min):
    #                         field_min= results[1] 
    #                         pdf_params[:,iy,ix] = parfit 
                            
    
    end_time = datetime.now()        
    log.info(" proc-ID "+str(process_id)+" ended at {}".format(end_time))
    print(" proc-ID "+str(process_id)+" ended at {}".format(end_time))
    print(" proc-ID "+str(process_id)+" elapsed time {}".format(end_time-start_time))
    return pdf_params

# ----------------------------------------------------------------------
#                ***  Write NetCDF file  ***
# Purpose: write the output parameters of the PDF fitting to a NetCDF file
# ----------------------------------------------------------------------
def write_netcdf(path,ntime,latitude,longitude, kappa1,loc1,lam1):

    nlon = len(longitude)
    nlat = len(latitude)
    time = ds0.variables["time"]
    
    fit_name =["skew", "weib","expweib", "laplace"]   
    month=[ "jan", "feb", "march", "april", 'may', 'june',"july", "aug", "sep","octo","nov","dec"] # If necessary for month wise analysis 
    exp_fit= target_pdf
    var_name= str.lower(input_var)+"_anom" 
    root_dir= "/home/mghani/atm_data/outputs_pdf/"   # change path/dir/folder  from the local computer 
    exp_folder = root_dir+"ensemble_01"
    
    if not os.path.isdir(exp_folder):
        os.mkdir(exp_folder)
    
    var_path= exp_folder+"{}".format(var_name) 
    if not os.path.isdir(var_path):
        os.mkdir(var_path)
    file_name = "pdfParm_{}_{}_{}_{}".format(var_name,exp_fit,in_yr,end_yr)+".nc"
    path_outdata= var_path +"/"+file_name
   
    print(path_outdata)
    ncfile = Dataset(path_outdata, 'w', format='NETCDF3_CLASSIC')
    
    #--- Open a new NetCDF file in write mode (specify format NETCDF3)
    # Create netCDF dimensions (lat and longitude)
    ncfile.createDimension("lat",nlat)
    ncfile.createDimension("lon",nlon)
    #ncfile.createDimension("time",ntime)
    
    # Create a netCDF variable (coordinate variable and variable).
    lat_nc  = ncfile.createVariable("latitude",np.float32,("lat",))
    lat_nc[:]  = latitude.values
    lat_nc.units  = latitude.attrs['units']

    lon_nc = ncfile.createVariable("longitude", np.float32,("lon",))
    lon_nc[:] = longitude.values
    lon_nc.units = longitude.attrs['units']
    
    #alpha = ncfile.createVariable("alpha", np.float32,("time","lat","lon" ))
    #print(np.shape(alpha[:]))  
    #alpha[:,:,:] = alpha1.values 
    #alpha.long_name = alpha1.attrs['long_name']

    kappa = ncfile.createVariable("kappa",np.float32,("lat","lon"))
    kappa[:,:] = kappa1.values
    kappa.long_name = kappa1.attrs['long_name']

    loc = ncfile.createVariable("loc",np.float32,("lat","lon"))
    loc[:,:] = loc1.values
    loc.long_name = loc1.attrs['long_name']

    lam = ncfile.createVariable("lambda",np.float32,("lat","lon"))
    lam[:,:] = lam1.values
    lam.long_name = lam1.attrs['long_name'] 
    
    ncfile.close()


#  Main function for reading input files and parallelization
# ----------------------------------------------------------------------

if __name__=='__main__':

    time_start = time.time()
    date_start = datetime.now()
    
    in_yr= sys.argv[1]
    end_yr =sys.argv[2]
    in_mon= sys.argv[3]
    out_mon=sys.argv[4]
    input_var=sys.argv[5] 
    target_pdf=sys.argv[6] #log=None  
    #input_var2=sys.argv[7]  #  input_var2 = var_vector2 for wind_vectors 
      
    start_day=date(int(in_yr),int(in_mon),1)
    start_date= start_day.strftime("%Y-%m-%d")
    end_day=date(int(end_yr),int(out_mon),31)
    print(start_date)
    end_date= end_day.strftime("%Y-%m-%d")
    print(end_date)

    #---> Introducing a log file for the program output ---
    file_log= "pdf_{}_{}_{}.log".format('%02s' % in_yr,'%02s' % end_yr,input_var) 
    logging.basicConfig(filename=file_log,level=logging.DEBUG,format = '%(levelname)s: %(message)s',filemode='w') 
    log= logging.getLogger(file_log)
    log.setLevel(logging.DEBUG)
    
    log.info(" ------------------------------")
    log.info("| PDF DISTRIBUTION FITTING")
    log.info("| python version="+sys.version)
    log.info("| Start time {}".format(date_start))
    log.info(" ------------------------------")
    log.info("")   
    #<--- Log file ends ---

    print(" ------------------------------")
    print(" | PDF distribution fit ")
    print(" | python version="+sys.version)
    print(" | Start Date:", date_start)
    print(" ------------------------------")
     
    yr_mon_list = []
    for r in range(int(in_yr),int(end_yr)+1):
        for s in range (int(in_mon),int(out_mon)+1):
             yr_mon_list.append('%02i' %r +'m' +'%02i' %s)
    print("yr-month=",yr_mon_list)
    
    #--- Reading Input files from input years range 

    files_dir= "/home/mghani/atm_data/ecmwf_ensemble/temp/" # Change with local target path/directory/ for input files
    
    files_list= sorted(list())
    for (path, dirs, files) in os.walk(files_dir):
        files_list+=[file for file in files]

    target_files= list()
    for f in files_list:        # crd_lat= latitude
        # crd_lon= longitude
        for g in range(len(yr_mon_list)):
            if f.startswith(("ECMWF_ENS_BULKTAIR_y"+yr_mon_list[g])):
               target_files.append(f)
    
    input_files= sorted(target_files)

    i=0
    files =[]
    for i in range(len(input_files)):
        fl= files_dir+ input_files[i]
        files.append(fl)


    input_array= xr.open_mfdataset(files,combine='nested', concat_dim=['time'], )
    input_array= input_array[input_var]
    # print inpout array shape
    print("input_array_shape=", np.shape(input_array))
    
    input_ar= np.mean(input_array, axis=0)
    print("input_ar_shape=", np.shape(input_ar))
    input_arr= xr.DataArray(input_ar, dims=['step','latitude','longitude'])
    print("input_array_shape=", np.shape(input_arr)) 

    print(len(files))
    print("1st_input_file=",files[0])
    print("Last_input_file=",files[-1])
    
    with xr.open_dataset(files[-1]) as ds_mesh:
        print(ds_mesh[input_var].shape)
        longitude = ds_mesh.variables["longitude"][:] 
        latitude  = ds_mesh["latitude"][:] 
        crd_lon =ds_mesh.variables["longitude"][:] 
        crd_lat = ds_mesh.variables["latitude"][:] 
        nlon = len(crd_lon)
        nlat = len(crd_lat)
        
    nday= len(input_ar['step'])
    crd_time=range(nday)
    print("nday=",nday)
   
    #---> Fixed conditions ---
    var_vector1=input_var
    #var_vector2=input_var2   # Use var_vector2 from input_var2, when compute wind_vectors (U,V) 
    lvar_vector=False          # lvar_vector= True, when we read wind vectors(U,V) inside worker_function 
    wind_vectors=False        #    wind_vectors=True, if we load and windv together to compute wind amplitude inside main function here
    pdf_fit= target_pdf       # target_pdf is retreived with the user input as pdf name 
 
    
    #-- processing for anomaly time series          
    def anom_tseries(xarr,start_yr_mon_day) :
        dt= pd.date_range(start=start_yr_mon_day, end=end_date, freq='6H') # freq=daily 
        ds_xarr= xr.DataArray(xarr,dims=['time','lat','lon'], coords=dict(time=dt))
        ds_cl =ds_xarr.groupby('time.month').mean('time') 
        ds_an= ds_xarr.groupby('time.month')-ds_cl
        return ds_an
    
# -----------------------------------------
    ##---> wind_vectors =True, and we use this section to read windu and windv together here 
    ## to compute wind amplitude, then broadcast values of input variable (as ds_anom variable) for the field values 
    ## at the bgeining of programme if this section is used <---|

    # windu= xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32),dims=['time','lat','lon'])
    # windv= xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32),dims=['time','lat','lon']) 
    # var_xarr=xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32), dims=['time','lat','lon']) 
   
    # iday=0
    # for iday in range(nday):
    #     file= xr.open_dataset(files[iday])
    #     file_list[iday] = fls[iday]
    #     var_xarr[iday,:,:]=  np.mean (file.variables[input_var][:,:,:], axis=0) #multiply with *0.01, when read MSL
    #     #windu[iday,:,:]= np.mean(file.variables["U10M"][:,:,:], axis=0) 
    #     #windv[iday,:,:]= np.mean(file.variables["V10M"][:,:,:], axis=0)  
    #     iday+=1
    # wndm = np.sqrt((windu[:,:,:]**2 + windv[:,:,:]**2))
  
    # if(wind_vectors):
    #     ds_anom1 = anom_tseries(windu, start_date)
    #     ds_anom2 = anom_tseries(windv, start_date)
    #     ds_anom  = np.sqrt((ds_anom1[:,:,:]**2 + ds_anom2[:,:,:]**2)) 
    # else:
    #     ds_anom = anom_tseries(var_xarr, start_date) #start date is always the first day of the time series

    pdf_var = input_arr[:,:,:] 
    pdf_len= len(pdf_var)
    print("var_shape=", np.shape(pdf_var))
    print("len_var", (pdf_len))
     
    #<---end--- wind_vectors (winud, winndv) loading and processing together 

    log.info(f" Configuration Paramiters:")
    log.info(f" Time Series:   year_start = {in_yr}")
    log.info(f'                year_end = {end_yr}')
    log.info(f'                month_start = {in_mon}')
    log.info(f'                month_end = {out_mon}')
    log.info(f'                input_variable = {input_var}')
    log.info(f'                fitcurve  = {pdf_fit}')
    
    #print(f" Configuration Paramiters:")
    print(f" Time Series:   year_start = {in_yr}")
    print(f'                year_end = {end_yr}')
    print(f'                month_start = {in_mon}')
    print(f'                month_end = {out_mon}')
    print( f'               input_variable = {input_var}')
    print(f' Fitting-Curve: fitcurve = {pdf_fit}')
    
   # -----------------------------------------
    # Parallelize the execution of a function distributing for input data/array  
    # -----------------------------------------
    #--> ncores=!!, we can change and check number of core required in our HPC system 
    avail_cores= cpu_count()
    # Use the available cores or a specified number, whichever is smaller
    ncores = min(avail_cores, 6)  # Adjust 6 to your desired number if needed

    print(f"Using {ncores} cores for processing")
    log.info(f"Using {ncores} cores for processing")
    log.info(f' Divide your spatial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
    #print(f' Divide your spacial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
    
    ilon = np.array_split(range(nlon), ncores)
    mylist =  [[ilon[i]] for i in range(ncores)]
    print(mylist[1][0][0])

    [mylist[i].append(crd_lon)  for i in range(ncores)]
    [mylist[i].append(crd_lat)  for i in range(ncores)]
    [mylist[i].append(crd_time)  for i in range(ncores)]
    #[mylist[i].append(files) for i in range(ncores)]
    #[mylist[i].append(var_vector1) for i in range(ncores)]
    
    log.info(f" Shape of the Pool data {np.shape(mylist)}")
    print(f" Shape of the Pool data (ncores x nelements) {np.shape(mylist)}")
    print(f" mylist[0][0] is index    Size = {len(mylist[0][0])}")
    print(f" mylist[0][1] is crd_lon Size = {len(mylist[0][1])}")
    print(f" mylist[0][2] is crd_lat      value = {mylist[0][2]}")
    #print(f" mylist[0][4] is crd_time      value = {mylist[0][3]}")
    #print(f" mylist[0][4] is var_vector  value = {mylist[0][5]}") 
   
    for icore in range(ncores):
        log.info(f" Assign {len(mylist[icore][0])}x{len(mylist[icore][1])} grid points to core{icore}")
        print(f" Assign {len(mylist[icore][0])}x{len(mylist[icore][1])} grid points to core{icore}")

    ############################################################################
    # p.map distributes function along the list
    ############################################################################
    
    with Pool(processes=ncores) as p:     
        results = p.map(worker_func, mylist)
        #print('Results (pool):\n', results)

    #---> Here, results are assigned to empty arrays 
    ds0 = xr.open_dataset(files[0])

    if (len(ds0.dims) == 3):

        kappa = xr.DataArray(np.zeros(((nlat,nlon)), dtype=np.float32), coords=[latitude, longitude], dims=['latitude','longitude'])
        #alpha = xr.DataArray(np.zeros(((nl
        loc   = xr.DataArray(np.zeros(((nlat,nlon)),dtype=np.float32), coords=[latitude, longitude], dims=['latitude','longitude'])
        lam   = xr.DataArray(np.zeros(((nlat,nlon)),dtype=np.float32 ), coords=[latitude, longitude], dims=['latitude','longitude'])

    else:
        print(f" Error! Dimension not correct")
   
    kappa.attrs['long_name'] = 'shape-kappa parameter'
    loc.attrs['long_name']   = 'location parameter'
    lam.attrs['long_name']   = 'scale-lambda parameter'
    
    out_path ="/path/dir"  # Change path/dir for save the output files 
    
    ilon=0
    if (len(ds0.dims) == 3):
        for result in results:
            #print('Result ',ilon,' shape(result):\n', np.shape(result))
            nlon_pool = len(result[0,0,:])
            #alpha[:,ilon:ilon+nlon_pool] = result[0,:,:]
            kappa[:,ilon:ilon+nlon_pool] = result[0,:,:]
            loc[:,ilon:ilon+nlon_pool]   = result[1,:,:]
            lam[:,ilon:ilon+nlon_pool]   = result[2,:,:]
            ilon += nlon_pool
        write_netcdf(out_path, nday,latitude,longitude,kappa,loc,lam) # We can add another paramter(alpha) depend on the PDF analysis  
 
    time_end = time.time()
    elapsed_time = (time_end - time_start)
    log.info(" elapsed time {}".format(elapsed_time)+"sec = {}".format(elapsed_time/60.)+"min")
    print(" elapsed time {}".format(elapsed_time)+ "sec = {}".format(elapsed_time/60.)+ " min")
    
    #--->log_path= change according to local path
    log_path="path/folder"  #  Select local path/folder to save the log file
    if not os.path.isdir(log_path):
        os.mkdir(log_path)
    src_logpath=os.getcwd()
    shutil.move(os.path.join(src_logpath,file_log),os.path.join(log_path, file_log))
