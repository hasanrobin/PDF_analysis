
from netCDF4 import Dataset
import numpy as np
import os

# ----------------------------------------------------------------------
#                  ***  write_fitParams_2D  ***
#
# Purpose:  Write Weibull fit parameters in a NETCDF file
# ----------------------------------------------------------------------
def write_fitParams_2D(path_outdata,file_outdata,longitude,latitude,
                       dim_lon,dim_lat,alpha1,kappa1,loc1,lam1):

	nlon = len(longitude)
	nlat = len(latitude)

	# Open a new NetCDF file in write mode (specify format NETCDF3)
	filepath_outdata = path_outdata+'/'+file_outdata
	#print(fileppath_outdata)
	ncfile = Dataset(filepath_outdata,'w', format='NETCDF3_CLASSIC')

	# Create netCDF dimensions (interv, lat and longitude)
	ncfile.createDimension(dim_lat,nlat)
	ncfile.createDimension(dim_lon,nlon)

	# Create a netCDF variable (coordinate variable and variable).
	lat_nc  = ncfile.createVariable(dim_lat,np.float32,(dim_lat,))
	lat_nc[:]  = latitude.values
	lat_nc.units  = latitude.attrs['units']

	lon_nc = ncfile.createVariable(dim_lon,np.float32,(dim_lon,))
	lon_nc[:] = longitude.values
	lon_nc.units = longitude.attrs['units']

	alpha = ncfile.createVariable("alpha",np.float32,(dim_lat,dim_lon))
	alpha[::] = alpha1.values
	alpha.long_name = alpha1.attrs['long_name']

	kappa = ncfile.createVariable("kappa",np.float32,(dim_lat,dim_lon))
	kappa[::] = kappa1.values
	kappa.long_name = kappa1.attrs['long_name']

	loc = ncfile.createVariable("loc",np.float32,(dim_lat,dim_lon))
	loc[::] = loc1.values
	loc.long_name = loc1.attrs['long_name']

	lam = ncfile.createVariable("lambda",np.float32,(dim_lat,dim_lon))
	lam[::] = lam1.values
	lam.long_name = lam1.attrs['long_name']

	# Close the file.
	ncfile.close()


# ----------------------------------------------------------------------
#                  ***  write_fitParams_2D  ***
#
# Purpose:  Write Weibull fit parameters in a NETCDF file
# ----------------------------------------------------------------------
def write_fitParams_3D(path_outdata,file_outdata,longitude,latitude,depth,
                       dim_lon,dim_lat,dim_depth,alpha1,kappa1,loc1,lam1):

	nlon   = len(longitude)
	nlat   = len(latitude)
	ndepth = len(depth)

	# Open a new NetCDF file in write mode (specify format NETCDF3)
	filepath_outdata = path_outdata+'/'+file_outdata
	#print(fileppath_outdata)
	ncfile = Dataset(filepath_outdata,'w', format='NETCDF3_CLASSIC')

	# Create netCDF dimensions (interv, lat and longitude)
	ncfile.createDimension(dim_lat,nlat)
	ncfile.createDimension(dim_lon,nlon)
	ncfile.createDimension(dim_depth,ndepth)

	# Create a netCDF variable (coordinate variable and variable).
	lat_nc  = ncfile.createVariable(dim_lat,np.float32,(dim_lat,))
	lat_nc[:]  = latitude.values
	lat_nc.units  = latitude.attrs['units']

	lon_nc = ncfile.createVariable(dim_lon,np.float32,(dim_lon,))
	lon_nc[:] = longitude.values
	lon_nc.units = longitude.attrs['units']

	depth_nc = ncfile.createVariable(dim_depth,np.float32,(dim_depth,))
	depth_nc[:] = depth.values
	depth_nc.units = depth.attrs['units']

	alpha = ncfile.createVariable("alpha",np.float32,(dim_depth,dim_lat,dim_lon))
	alpha[::] = alpha1.values
	alpha.long_name = alpha1.attrs['long_name']

	kappa = ncfile.createVariable("kappa",np.float32,(dim_depth,dim_lat,dim_lon))
	kappa[::] = kappa1.values
	kappa.long_name = kappa1.attrs['long_name']

	loc = ncfile.createVariable("loc",np.float32,(dim_depth,dim_lat,dim_lon))
	loc[::] = loc1.values
	loc.long_name = loc1.attrs['long_name']

	lam = ncfile.createVariable("lambda",np.float32,(dim_depth,dim_lat,dim_lon))
	lam[::] = lam1.values
	lam.long_name = lam1.attrs['long_name']

	# Close the file.
	ncfile.close()











# ----------------------------------------------------------------------
#                  ***  write_moments_2D  ***
#
# Purpose:   Write moments in a NETCDF file Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis
# ----------------------------------------------------------------------
def write_moments_2D(path_outdata,file_outdata,longitude,latitude,dim_lon, dim_lat
                        ,mean1,variance1,skew1,kurtosis1):

	nlon = len(longitude)
	nlat = len(latitude)

	# Open a new NetCDF file in write mode (specify format NETCDF3)
	filepath_outdata = path_outdata+'/'+file_outdata
	# print(filepath_outdata)
	ncfile = Dataset(filepath_outdata,'w', format='NETCDF3_CLASSIC')

	# Create netCDF dimensions (interv, lat and longitude)
	ncfile.createDimension(dim_lat,nlat)
	ncfile.createDimension(dim_lon,nlon)

	# Create a netCDF variable (coordinate variable and variable).
	lat_nc  = ncfile.createVariable(dim_lat,np.float32,(dim_lat,))
	lat_nc[:]  = latitude.values
	#lat_nc.units  = latitude.attrs['units']

	lon_nc = ncfile.createVariable(dim_lon,np.float32,(dim_lon,))
	lon_nc[:] = longitude.values
	#lon_nc.units = longitude.attrs['units']

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
#                  ***  write_moments_3D  ***
#
# Purpose:  Write moments in a NETCDF file Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis
# ----------------------------------------------------------------------
def write_moments_3D(path_outdata,file_outdata,longitude,latitude,depth,
                     dim_lon,dim_lat,dim_depth,mean1,variance1,skew1,kurtosis1):

	nlon   = len(longitude)
	nlat   = len(latitude)
	ndepth = len(depth)

	# Open a new NetCDF file in write mode (specify format NETCDF3)
	filepath_outdata = path_outdata+'/'+file_outdata
	#print(fileppath_outdata)
	ncfile = Dataset(filepath_outdata,'w', format='NETCDF3_CLASSIC')

	# Create netCDF dimensions (interv, lat and longitude)
	ncfile.createDimension(dim_lat,nlat)
	ncfile.createDimension(dim_lon,nlon)
	ncfile.createDimension(dim_depth,ndepth)

	# Create a netCDF variable (coordinate variable and variable).
	lat_nc  = ncfile.createVariable(dim_lat,np.float32,(dim_lat,))
	lat_nc[:]  = latitude.values
	lat_nc.units  = latitude.attrs['units']

	lon_nc = ncfile.createVariable(dim_lon,np.float32,(dim_lon,))
	lon_nc[:] = longitude.values
	lon_nc.units = longitude.attrs['units']

	depth_nc = ncfile.createVariable(dim_depth,np.float32,(dim_depth,))
	depth_nc[:] = depth.values
	depth_nc.units = depth.attrs['units']

	mean = ncfile.createVariable("mean",np.float32,(dim_depth,dim_lat,dim_lon))
	mean[::] = mean1.values
	mean.long_name = mean1.attrs['long_name']

	variance = ncfile.createVariable("variance",np.float32,(dim_depth,dim_lat,dim_lon))
	variance[::] = variance1.values
	variance.long_name = variance1.attrs['long_name']

	skew = ncfile.createVariable("skew",np.float32,(dim_depth,dim_lat,dim_lon))
	skew[::] = skew1.values
	skew.long_name = skew1.attrs['long_name']

	kurtosis = ncfile.createVariable("kurtosis",np.float32,(dim_depth,dim_lat,dim_lon))
	kurtosis[::] = kurtosis1.values
	kurtosis.long_name = kurtosis1.attrs['long_name']

	# Close the file.
	ncfile.close()
