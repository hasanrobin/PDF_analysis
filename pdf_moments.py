#!/usr/bin/python3
# Python package requirements

import os
import sys
import time
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, date, timedelta
import xarray as xr
import getopt
import glob 
import scipy
from scipy import stats
from multiprocessing import Pool, cpu_count


def worker_func(mylist):

    process_id = os.getpid()
    start_time = datetime.now()
    #log.info(" proc-ID "+str(process_id)+" started at {}".format(start_time))
    print(" proc-ID "+str(process_id)+" started at {}".format(start_time))

    nlon_sub = len(mylist[0])
    crd_lon      = mylist[1]
    crd_lat      = mylist[2]

    # print(crd_lon)
    # print(crd_lat) 
    # print(nlon_sub)
    ds0 = xr.open_dataset(input_pdf, engine="netcdf4")
    # print(ds0)
    nlat     = len(crd_lat)
    nlon     = len(crd_lon)
    
    if (len(ds0.dims) == 2):
        ndepth = 1
    
    lambd = xr.DataArray(np.zeros((nlat,nlon_sub)))
    kappa = xr.DataArray(np.zeros((nlat,nlon_sub)))
    alpha = xr.DataArray(np.zeros((nlat,nlon_sub)))
    locat = xr.DataArray(np.zeros((nlat,nlon_sub)))

    moments = xr.DataArray(np.zeros((4,nlat,nlon_sub)))
    #moments[:,:,:] = np.nan
    #print(mylist[0][0],mylist[0][nlon_sub-1]+1)

    if (len(ds0.dims) == 2):
        lambd[:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['lambda']
        kappa[:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['kappa']
        #alpha[:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['alpha']
        locat[:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['loc']
    elif (len(ds0.dims) == 3):
        lambd[:,:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['lambda']
        kappa[:,:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['kappa']
        alpha[:,:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['alpha']
        locat[:,:,:] = ds0.isel(lon=slice(mylist[0][0],mylist[0][nlon_sub-1]+1)).load()['loc']
    print(lambd.sizes)
    # print(kappa.sizes)
    # print(alpha.sizes)
    # print(locat.sizes)
    
    for iz in range(ndepth):
        for iy in range(nlat):
            for ix in range(nlon_sub):
                #if not isnan(kappa[iy,ix]):             
                    #moments[0,iy,ix], moments[1,iy,ix], moments[2,iy,ix], moments[3,iy,ix] = scipy.stats.weibull_min.stats(kappa[iy,ix],loc=0,scale=lambd[iy,ix], moments='mvsk')
                     #moments[0,iy,ix], moments[1,iy,ix], moments[2,iy,ix], moments[3,iy,ix] = scipy.stats.laplace_asymmetric.stats(kappa[iy,ix],loc=locat[iy,ix],scale=lambd[iy,ix], moments='mvsk')
                     moments[0,iy,ix], moments[1,iy,ix], moments[2,iy,ix], moments[3,iy,ix] = scipy.stats.skewnorm.stats(kappa[iy,ix],loc=locat[iy,ix],scale=lambd[iy,ix], moments='mvsk')

    end_time = datetime.now()
    diff_time = (end_time - start_time).seconds
    # log.info(" proc-ID "+str(process_id)+" completed at {}".format(end_time)+" after {}".format(diff_time)+"sec = {}".format(diff_time/60.)+"min")
    print(" proc-ID "+str(process_id)+" completed at {}".format(end_time)+" after {}".format(diff_time)+"sec = {}".format(diff_time/60.)+"min")

    return moments

# ----------------------------------------------------------------------
#                  ***  write_moments_2D  ***
# Purpose:   Write moments in a NETCDF file Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis
# ----------------------------------------------------------------------
def write_moments_2D(path_outdata,file_outdata,longitude,latitude,
                       mean1,variance1,skew1,kurtosis1):

    nlon = len(longitude)
    nlat = len(latitude) 
    dim_lat = "lat"
    dim_lon= "lon" 
	# Open a new NetCDF file in write mode (specify format NETCDF3)
    filepath_outdata = path_outdata+'/'+file_outdata
	# print(filepath_outdata)
    ncfile = Dataset(filepath_outdata,'w', format='NETCDF3_CLASSIC')

	# Create netCDF dimensions (interv, lat and longitude)
    ncfile.createDimension("lat",nlat)
    ncfile.createDimension("lon",nlon)

	# Create a netCDF variable (coordinate variable and variable).
    lat_nc  = ncfile.createVariable("laitude",np.float32,("lat",))
    lat_nc[:]  = latitude.values
    lat_nc.units  = latitude.attrs['units']

    lon_nc = ncfile.createVariable("longitude",np.float32,("lon",))
    lon_nc[:] = longitude.values
    lon_nc.units = longitude.attrs['units']

    mean = ncfile.createVariable("mean",np.float32,(dim_lat,dim_lon))
    mean[::] = mean1.values
    mean.long_name = mean1.attrs['long_name']

    variance = ncfile.createVariable("variance",np.float32,(dim_lat,dim_lon))
    variance[::] = variance1.values
    variance.long_name = variance1.attrs['long_name']

    skew = ncfile.createVariable("skew",np.float32,(dim_lat,dim_lon))
    skew[::] = skew1.values
    skew.long_name = skew1.attrs['long_name']

    kurtosis = ncfile.createVariable("kurtosis",np.float32,(dim_lat,dim_lon))
    kurtosis[::] = kurtosis1.values
    kurtosis.long_name = kurtosis1.attrs['long_name']

	# Close the file.
    ncfile.close()


# ----------------------------------------------------------------------
if __name__=='__main__':

    time_start = time.time()
    date_start = datetime.now()
    
    print(" ------------------------------")
    print(" | MOMENTS OF PDF fitting DISTRIBUTION PARAMETERS")
    print(" | python version="+sys.version)
    print(" | Start Date:", date_start)
    print(" ------------------------------")
    print("")
    #log=None
     
    input_var= sys.argv[1]

    # Read script args
    var_name = None
    try:
        #options, rem = getopt.getopt(sys.argv[1:], 'd:n:e:h', ['idexp=','ncores=','path_exp=', 'help'])
        options, rem = getopt.getopt(sys.argv[1:], 'v:', ['var_name=' 'help'])
        print(options,rem)
        for opt,arg in options:
            if opt in ('-v', '--var_name'):
                var_name = arg
            elif opt in ('-h', '--help'):
                #printHelp(log)
                sys.exit(1)

    except getopt.GetoptError:
        print("ERROR! no variable name!")
        #printHelp(log)
        sys.exit(1)

    print(f" Script Args: variable  = {var_name}")
    #print(f'              NCORES = {ncores}')
    #print(f'              path_exp = {path_exp}') 

    var_name= str.lower(var_name)+"_anom"
    fitcurve =["skew","weib", "laplace"] 
    
    print(f" Fixed Parameters:")
    print(f'                var_vector1 = {var_name}')
    print(f' Fitting-Curve: fitcurve = {fitcurve[0]}')
    
    #-- Read variable latitude and longitude arrays
    files_list=glob.glob("path/dir/"+"*.nc")   # Select the path/dir for an input file
    with xr.open_dataset(files_list[0]) as ds0:
        #print(ds0)
        longitude= ds0.variables["lon"]
        latitude= ds0.variables["lat"]
        crd_lon= longitude
        crd_lat = latitude 
        nlat=  len(latitude)
        nlon = len(longitude)

        if (len(ds0.dims) == 2):
           ndepth = 1
        
    
    #If the experiment out folder does not exist, create it
  
    exp_out_dir= '/path/dir/'   # Change the path/dir to local dir
    exp_var_path= exp_out_dir+"{}".format(var_name)
    print(exp_var_path)
    if not os.path.isdir(exp_var_path):
        print("ERROR! PATH "+ exp_var_path+" does not exist!")
        os.mkdir(exp_var_path)
        sys.exit(0)

    path_fitdata = exp_out_dir+ '/'+var_name
    input_pdf = exp_var_path +"/"+"pdfParm_{}_{}_2006_2020.nc".format(var_name, fitcurve[0] )
    file_pdfMom = "fitmom_{}_{}_06_20.nc".format(var_name, fitcurve[0])
    print(input_pdf) 

    # -----------------------------------------
    # Parallelize the execution of a function distributing the input data across processes
    # -----------------------------------------
    # log.info(f' Divide your spacial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
    
    ncores=36
    print(f' Divide the spatial domanin {len(longitude)} x {len(latitude)} for {ncores} cores')
   
    ilon = np.array_split(range(nlon), ncores)
    # print(ilon)
    mylist = [[ilon[i]] for i in range(ncores)]

    #[mylist[i].append(path_fitdata) for i in range(ncores)]
    [mylist[i].append(crd_lon)      for i in range(ncores)]
    [mylist[i].append(crd_lat)      for i in range(ncores)]
    
    # Start ncores worker processes and do the computation
    with Pool(processes=ncores) as p:
        results = p.map(worker_func, mylist)
        #print('Results (pool):\n', results)

    # Now the results are assigned to empty arrays
    ds0 = xr.open_dataset(input_pdf, engine="netcdf4")
    if (len(ds0.dims) == 2):
        mean     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        variance = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        skew     = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
        kurtosis = xr.DataArray(np.zeros((nlat,nlon)), coords=[latitude, longitude], dims=['latitude','longitude'])
    
    else:
        print(f" Error! Dimension not correct")

    mean.attrs['long_name']     = 'mean'
    variance.attrs['long_name'] = 'variance'
    skew.attrs['long_name']     = 'skew'
    kurtosis.attrs['long_name'] = 'kurtosis'

    ilon=0
    if (len(ds0.dims) == 2):
        for result in results:
            # print('Result ',ilon,' shape(result):\n', np.shape(result))
            nlon_pool = len(result[0,0,:])
            # print(nlon_pool)
            mean[:,ilon:ilon+nlon_pool]     = result[0,:,:]
            variance[:,ilon:ilon+nlon_pool] = result[1,:,:]
            skew[:,ilon:ilon+nlon_pool]     = result[2,:,:]
            kurtosis[:,ilon:ilon+nlon_pool] = result[3,:,:]
            ilon += nlon_pool

        write_moments_2D(exp_var_path,file_pdfMom,longitude,latitude,
                            mean,variance,skew, kurtosis)


    time_end = time.time()
    elapsed_time = (time_end - time_start)
    print(exp_var_path+'/'+file_pdfMom)
    print(" Elapsed time {}".format(elapsed_time)+"sec = {}".format(elapsed_time/60.)+"min")

