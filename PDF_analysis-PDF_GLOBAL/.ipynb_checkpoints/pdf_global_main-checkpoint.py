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
import warnings 
import  pandas as pd
sys.path.insert(0,"/work/cmcc/mg13420/plot_exercises/jup_nbooks/method_function/")
from file_directory import* 
os.system("hostname")  
#import statsmodels.tsa.seasonal as tsa 
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
    #print(crd_time)
    ndays    =  len(crd_time)
    mlat     =  len(crd_lat)
    mlon     =  len(crd_lon)
    
# -----------------------------------------------------
# Define/upload the variable 'field'/values to be fitted
# ----------------------------------------------------- 
    
    if 'latitude' in ds0.variables:
        field = xr.DataArray(np.zeros(((ndays, mlat,nlon_sub)), dtype=np.float32),dims=['time', 'latitude','longitude'])
        field[:,:,:] = anom_var.isel(longitude=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[:,:,:]
    if 'lat' in ds0.variables:
        field = xr.DataArray(np.zeros(((ndays, mlat,nlon_sub)), dtype=np.float32),dims=['time', 'lat','lon']) 
        field[:,:,:] = anom_var.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[:,:,:] # Use this if we broadcast input values from the main function section
    
    field =np.nan_to_num(field, nan=-500)
    #print("field_mean", (field[0].values) )
    
    # add the statsmodels library to compute the residuals using parallel processing 
    # resid_arr = xr.DataArray(np.zeros(((ndays, mlat,nlon_sub)), dtype=np.float64),dims=['time', 'lat','lon'])
    # # for iday in range(ndays):
    # for ilat in range(mlat):
    #     for ilon in range(nlon_sub):
    #         if not np.isnan(field[0,ilat,ilon]):
    #             #continue
    #             result = tsa.seasonal_decompose(field[:,ilat,ilon], model='additive', period= 7 )

    #             resid =result.resid 
    #             resid = np.nan_to_num(resid, nan=np.nanmean(resid))
    #             resid_arr[:,ilat,ilon] = result.resid

    # field1 = xr.DataArray(resid_arr, dims=['time','lat','lon'], )
    # print("Values", np.nanmean(resid_arr[:]) )

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
    #             ds = xr.open_dataset(files[0])
    #             field[iday,:,:] = ds.isel(longitude=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()[var_vector1][0]# multiply with *0.01 for MSLP
    #             #print(f"Shape of field[{iday},:,:]: {np.shape(field[iday,:,:])}")
    #            iday+=1
            
    
# Define the FIT parameters, PDF  statistics
# -----------------------------------------
    
    pdf_params = xr.DataArray(np.zeros((nparams, mlat,nlon_sub), dtype=np.float32), dims=['param', 'latitude','longitude'])
    pdf_params[:,:,:] = np.nan
    #print(np.shape(pdf_params))

    #  fit the PDF to the data depending on the number of parameters 
    if (nparams == 2):
        #for i in range(nday):
        for iy in range(mlat):
            for ix in range(nlon_sub):
                #if not isnan(field[0,iy,ix]):
                if (field[0,iy,ix]==-500):
                    continue
                pdf_params[:,iy,ix] = list(pdf_name.fit(field[:,iy,ix]))
                
    # Provide initial guesses for the parameters
    initial_guesses = {'skewnorm': (0.1,np.nanmean(field), np.nanstd(field))}
    if (nparams == 3):
        for iy in range(mlat):
            for ix in range(nlon_sub):
                #if not isnan(field[0,iy,ix]):
                if (field[0,iy,ix]==-500):
                    continue
                pdf_params[:,iy,ix] = list(pdf_name.fit(field[:,iy,ix],0.005))  # 0.01 works fine for the initial guess
        
    # elif (nparams ==3):
    # # # #---This skewnormal PDF function is adapted with the Nelder-med algorithm, if we have convergent issue in shape parameter distribution  
    #     if (pdf_fit == "skewnorm"):
    #         for iy in range(mlat):
    #             for ix in range(nlon_sub):
    #                 if (field[0,iy,ix]==-500):
    #                     continue
    #                 #if np.isfinite(field[0,iy,ix]):
    #                 n=3
    #                 # We can try with different number set 
    #                 par_min=  [-15,field.min(),0.1] 	
    #                 par_max=  [15, field.max(), 10] 
    #                 rand_parm =np.random.uniform(low=par_min, high =par_max, size =(n,3))
    #                 field_min=10**10
    #                 for ipar in range(len(rand_parm)):
    #                     results = list(pdf_name.fit(field[:,iy,ix],rand_parm[ipar][0],loc=rand_parm[ipar][1],scale=rand_parm[ipar][2] ))
    #                     parfit = results[0]
    #                     if (results[1]< field_min):
    #                         field_min= results[1] 
    #                         pdf_params[:,iy,ix] = parfit
    #                     ipar+=1
    #     # else:
        #     for iy in range(mlat):
        #         for ix in range(nlon_sub):
        #             if np.isfinite(field[0,iy,ix]):
        #                 pdf_params[:,iy,ix] = list(pdf_name.fit(field[:,iy,ix], )) 

    end_time = datetime.now()        
    log.info(" proc-ID "+str(process_id)+" ended at {}".format(end_time))
    print(" proc-ID "+str(process_id)+" ended at {}".format(end_time))
    print(" proc-ID "+str(process_id)+" elapsed time {}".format(end_time-start_time))
    return pdf_params

# ----------------------------------------------------------------------
#                ***  Write NetCDF file  ***
# Purpose: write the output parameters of the PDF fitting to a NetCDF file
# ----------------------------------------------------------------------
def write_netcdf(path,ntime,latitude,longitude,kappa1,loc1,lam1):

    nlon = len(longitude)
    nlat = len(latitude)
    time = ds0.variables["time"]
    
    fit_name =["skew", "weib","expweib", "laplace"]   
    month=[ "jan", "feb", "march", "april", 'may', 'june',"july", "aug", "sep","octo","nov","dec"] # If necessary for month wise analysis 
    exp_fit= target_pdf
    var_name= str.lower(input_var) 
    root_dir= "/work/cmcc/mg13420/plot_exercises/"   # change path/dir/folder  from the local computer 
    exp_folder = root_dir+"sst_pdf-exps/"  # Change folder name according to the experiment
    
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

    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        warnings.filterwarnings("ignore")

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
     
    #--- Reading Input files from input years range 
    # files_dir= "/data/cmcc/mg13420/model_atmos/fc_ensemble/temp/" # Change with local target path/directory/ for input files
    #files_dir="/data/cmcc/mg13420/sst_data/model/"
    files_dir="/data/cmcc/mg13420/sst_data/sst_IBI_0.05deg/"
    #files_dir="/data/cmcc/mg13420/sst_data/obs_sst_IBI_0.05deg/"
    
    atm_comb=  target_years_files(files_dir,in_yr, end_yr, in_mon, out_mon, target_files='sst', era5_file=None)
        
    #atm_comb=  input_yr_files(files_dir,in_yr, end_yr, 'sst', None)
    print("var_data", len(atm_comb['time']))

    with xr.open_dataset(glob(files_dir+"*.nc")[0]) as ds0:
        print("single file", ds0[input_var].shape)
        latitude=  atm_comb.variables["lat"][:]
        longitude=  atm_comb.variables["lon"][:]
        # longitude = ds0["longitude"][:] 
        # latitude  = ds0["latitude"][:] 
        # crd_lon =ds0.variables["longitude"][:] 
        # crd_lat = ds0.variables["latitude"][:] 
        crd_lon =ds0.variables["lon"][:] 
        crd_lat = ds0.variables["lat"][:] 
        
        nlon = len(crd_lon)
        nlat = len(crd_lat)

       
        if 'time' in ds0.variables:
            time_step = ds0.variables["time"]
        if 'step' in ds0.variables:
            time_step = ds0.variables["step"][:]

        if 'depth' in ds0.variables:
            depth = ds0.variables["depth"][0]
            print("Depth level=",depth.values) 
        else:
            depth = None
        
        nday= len(atm_comb['time'])*len(time_step)
        nday1= len(atm_comb['time'][:])
        crd_time=range(nday1)      
        print("Length of target time series=",nday1)
        print(nday, nlat, nlon)


    # if 'time' in ds0.variables:
    #     concat_dim = 'time'
    #     comb_files = xr.open_mfdataset(files, combine='nested', parallel=True, concat_dim=[concat_dim])
    #     print("Using time dimension for  concatenation")
    #     print("Combined files shape=", comb_files[input_var].shape)
    # if 'step' in ds_mesh.variables: 
    #     concat_dim = 'step'
    #     comb_files = xr.open_mfdataset(files, combine='nested', concat_dim=[concat_dim]) 
    #     print("Using step dimension for concatenation")
    #     print("Combined files shape=", comb_files[input_var].shape)
    
    # #comb_files=comb_files.swap_dims({'step':'time'})

   
    # #  Load targt variable and check if depth is present in the variable 

    # input_array = comb_files[input_var]

    # if 'depth' in input_array.dims:
    # # if 'depth' in comb_files.dims:
    #     # select first depth level
    #     if 'time' in input_array.dims:
    #         input_array =xr.DataArray(input_array[:,0,:,:], dims=['time','latitude','longitude'],  )
    #         print ("Depth free shape=", (input_array)) 
            
    #     else:
    #         input_array = xr.DataArray(input_array, dims=['step','latitude','longitude'])
        
    # else :
    #     if 'time' in input_array.dims:
    #         input_array = xr.DataArray(input_array, dims=['time','latitude','longitude'])
    #         print ("Darray (time, lat, lon) shape=", (input_array)) 
    #     else:
    #         input_array = xr.DataArray(input_array, dims=['step','latitude','longitude'])
    #         print ("Dataarray (step,lat,lon) shape=", np.shape(input_array)) 


    #--- Passing input variable to the worker function 
    #sst= atm_comb['analysed_sst'][:,:,:]  
    sst= atm_comb['thetao'][:,0,:,:]
    #sst = xr.DataArray(sst, dims=['time','latitude','longitude'], coords=[atm_comb['time'][:8], latitude, longitude]) 
    sst = xr.DataArray(sst, dims=['time','lat','lon'], coords=[atm_comb['time'], latitude, longitude]) 
    
    #lsm =xr.DataArray(lsm, dims=['time','lat','lon'], coords=[dt2, latitude, longitude]) 
    # FillValue = -54.53
    # pdf_var = np.where(sst==FillValue, np.nan, sst)
    pdf_var = sst.load()+ 273.15
    
    #---> Fixed conditions ---
    var_vector1=input_var
    #var_vector2=input_var2   # Use var_vector2 from input_var2, when compute wind_vectors (U,V) 
    lvar_vector=False          # lvar_vector= True, when we read wind vectors(U,V) inside worker_function 
    wind_vectors=False        #    wind_vectors=True, if we load and windv together to compute wind amplitude inside main function here
    pdf_fit= target_pdf       # target_pdf is retreived with the user input as pdf name 
 

    #--To control the number of parameters of the target PDF function
    pdf_fit= target_pdf          # target_pdf is retreived with the user input as pdf name 
    pdf_name= getattr(stats,pdf_fit)
    print("pdf_short_name=", pdf_name)
    rand_data= np.random.normal(size=100)   # random number for testing the PDF fitting)
   
    # check the number of parameters of the PDF
    pfit = pdf_name.fit(rand_data,)
    nparams = len(pfit[:3]) 
    pdf_fullname= pdf_name.__class__.__name__[:-4]
    print(f"PDF name : {pdf_fullname}  and  Number of PDF parameters= {nparams} ")
    log.info(f"PDF name :{pdf_fullname} and  Number of PDF parameters= {nparams} ") 
    
    #-- processing for anomaly time series          
    def anom_tseries(xarr,start_date, end_date):
        dt= pd.date_range(start=start_date, end=end_date, freq='D') # freq=daily 
        ds_xarr= xr.DataArray(xarr,dims=['time','lat','lon'], coords=dict(time=dt))
        ds_cl =ds_xarr.groupby('time.dayofyear').mean('time') 
        ds_an= ds_xarr.groupby('time.dayofyear')-ds_cl
        return ds_an
    
    anom_var = anom_tseries(pdf_var, start_date, end_date)
    print("Anomaly shape=", anom_var[0].values)
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

    # pdf_var = input_arr[:,:,:] 
    # pdf_len= len(pdf_var)
    # print("var_shape=", np.shape(pdf_var))
    # print("len_var", (pdf_len))
     
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
    ncores=  cpu_count()
    # Use the available cores or a specified number, whichever is smaller
    #ncores = min(avail_cores, 72)  # Adjust 6 to your desired number if needed

    print(f"Using {ncores} cores for processing")
    log.info(f"Using {ncores} cores for processing")
    log.info(f' Divide your spatial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
    print(f' Divide your spacial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
    
    ilon = np.array_split(range(nlon), ncores)
    mylist =  [[ilon[i]] for i in range(ncores)]
    mylist_np = np.array(mylist, dtype=object)
    [mylist[i].append(crd_lon)  for i in range(ncores)]
    [mylist[i].append(crd_lat)  for i in range(ncores)]
    [mylist[i].append(crd_time)  for i in range(ncores)]
    #[mylist[i].append(files) for i in range(ncores)]
    #[mylist[i].append(var_vector1) for i in range(ncores)]
    
    #log.info(f" Shape of the Pool data {np.shape(mylist)}")
    #print(f" Shape of the Pool data (ncores x nelements) {np.shape(mylist)}")
    print(f" mylist[0][0] is index    Size = {len(mylist[0][0])}")
    print(f" mylist[0][1] is crd_lon Size =  {len(mylist[0][1])}")
    # print(f" mylist[0][2] is crd_lat    value = {mylist[0][2]}")
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

    if (len(ds0.dims) == 2):
       
        loc   = xr.DataArray(np.zeros(((nlat,nlon)),dtype=np.float32),  dims=['latitude','longitude'])
        lam   = xr.DataArray(np.zeros(((nlat,nlon)),dtype=np.float32 ),  dims=['latitude','longitude'])

    elif(len(ds0.dims) >= 3):
        kappa = xr.DataArray(np.zeros(((nlat,nlon)), dtype=np.float32),  dims=['latitude','longitude']) 
        #alpha = xr.DataArray(np.zeros(((nday,nlat,nlon)), dtype=np.float32),  dims=['time','latitude','longitude'])
        loc   = xr.DataArray(np.zeros(((nlat,nlon)), dtype=np.float32),  dims=['latitude','longitude']) 
        lam   = xr.DataArray(np.zeros(((nlat,nlon)), dtype=np.float32),  dims=['latitude','longitude'])
    else:
        print(f" Error! Dimension not correct")
   
    kappa.attrs['long_name'] = 'shape-kappa parameter'
    loc.attrs['long_name']   = 'location parameter'
    lam.attrs['long_name']   = 'scale-lambda parameter'
    out_path ="/work/cmcc/mg13420/plot_exercises/"  # Change path/dir for save the output files 
    
    ilon=0
    if (len(ds0.dims) >= 3):
        for result in results:
            #print('Result ',ilon,' shape(result):\n', np.shape(result))
            nlon_pool =result.shape[2]
            nparams = result.shape[0] 
            if nparams == 2:
                kappa[:,ilon:ilon+nlon_pool] = result[0,:,:]
                loc[:,ilon:ilon+nlon_pool]   = result[1,:,:]

            elif nparams == 3:
                kappa[:,ilon:ilon+nlon_pool] = result[0,:,:]
                loc[:,ilon:ilon+nlon_pool]   = result[1,:,:]
                lam[:,ilon:ilon+nlon_pool]   = result[2,:,:]
            ilon += nlon_pool
            
        write_netcdf(out_path, nday,latitude,longitude,kappa,loc,lam) # We can add another paramter(alpha) depend on the PDF analysis  
 
    time_end = time.time()
    elapsed_time = (time_end - time_start)
    log.info(" elapsed time {}".format(elapsed_time)+"sec = {}".format(elapsed_time/60.)+" min")
    print(" elapsed time {}".format(elapsed_time)+ "sec = {}".format(elapsed_time/60.)+ " min")
    
    #--->log_path= change according to local path
    log_path=out_path  #  Select local path/folder to save the log file
    if not os.path.isdir(log_path):
        os.mkdir(log_path)
    src_logpath=os.getcwd()
    shutil.move(os.path.join(src_logpath,file_log),os.path.join(log_path, file_log))
