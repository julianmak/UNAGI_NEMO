#!/usr/bin/env python3

# JM: 25 Feb 2019
# script to generate the depth integrated epe
# requires already processed files of the isopyc_tave

import numpy as np
import glob, netCDF4, copy, sys

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = 
    """
    Script to generate the depth integrated epe, where by
    
      epe(x, y, sigma) = (g / 2 rho0) < z'z' > = (g / 2 rho0) (< z^2 > - < z >^2)
      
    where z = z(x, y, sigma) is the diagnosed isopycnal depths in density 
    (not buoyancy) co-ordinates
    """)
    
# fixed arguments
parser.add_argument("data_dir", type = str,
                    help = "specify data directory")
parser.add_argument("keys",     type = str, 
                    help = "specify string in filename to grab (e.g. \"*isodep_tave*\", note the quote marks, don't include the .nc extension)")

# optional arguments
parser.add_argument("--file_out", type = str, 
                    help = "specify output name (e.g. epe_zint_tave.nc)",
                    default = "epe_zint_tave.nc")
parser.add_argument("--lquery",
                    help = "print out the variables available",
                    action = "store_true")
parser.add_argument("--sigma_var", type = str,
                    help = "variable name for the density bin (default = gamma_a)",
                    default = "gamma_a")
parser.add_argument("--isodep_var", type = str,
                    help = "variable name for isodep (default = isodep_tave)",
                    default = "isodep_tave")
parser.add_argument("--isodepsq_var", type = str,
                    help = "variable name for isodepsq (default = isodepsq_tave)",
                    default = "isodepsq_tave")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# initialisation

# grab file names  
filenames = args.data_dir + args.keys + ".nc"
file_list = []
for file in glob.glob(filenames):
  file_list.append(file)
    
if len(file_list) == 0:
  sys.exit("no files grabbed in %s, are you in the right folder?" % filenames)

# sort it according to the timestamps
file_list.sort()

#--------------------------------------------------------
# cycle through the file lists and compile the moc
print("%g files found, cycling through them..." % len(file_list))

for i in range(len(file_list)):
  
  file = file_list[i].replace(args.data_dir, "") # strip out the data_dir
  print(" ")
  print("working in file = %g / %g" % (i + 1, len(file_list)))
  data = netCDF4.Dataset(file_list[i])

  if args.lquery:
    for name, variable in data.variables.items():
      for attrname in variable.ncattrs():
        if attrname == "standard_name":
          print("{} -- {}".format(name, getattr(variable, attrname)))

    data.close()

    print(" ")

    sys.exit("finished query, exiting gen_epe_tave...")
    
  else:
  
    # do extra averaging if need be
    
    if i == 0:
      nav_lon    = data.variables["nav_lon"][:, :]
      nav_lat    = data.variables["nav_lat"][:, :]
      sigma      = data.variables[args.sigma_var][:]
      isodep_tave   = data.variables[args.isodep_var][:, :, :]
      isodepsq_tave = data.variables[args.isodepsq_var][:, :, :]
      
      # only generate weights where there is a non-zero record
      wdep_tave     = np.where(isodep_tave != 0.0, 1, 0)
      isodep_tave   = data.variables[args.isodep_var][:, :, :]
      isodepsq_tave = data.variables[args.isodepsq_var][:, :, :]
    else:
      wdep_tave     += np.where(isodep_tave != 0.0, 1, 0)
      isodep_tave   += data.variables[args.isodep_var][:, :, :]
      isodepsq_tave += data.variables[args.isodepsq_var][:, :, :]

  data.close()

# average by the number of weights
isodep_tave   = np.where(wdep_tave > 0, isodep_tave   / wdep_tave, 0.0)
isodepsq_tave = np.where(wdep_tave > 0, isodepsq_tave / wdep_tave, 0.0)

# generate the epe field given by in density co-ordinates and vertically integrate
g, rho0 = 9.8, 1026.0

epe = 0.5 * (g / rho0) * (isodepsq_tave - isodep_tave * isodep_tave)

epe_zint = np.trapz(epe, sigma, axis = 0) # vertically integrate

#--------------------------------------------------------
# write the file

print("writing...")

# open a new netCDF file for writing.
ncfile = netCDF4.Dataset(args.data_dir + args.file_out, "w", format = "NETCDF4") 
ncfile.title = "(specific) epe file from UNAGI_R010 data"

# create the dimensions.
ncfile.createDimension("x", nav_lon.shape[1])
ncfile.createDimension("y", nav_lon.shape[0])
ncfile.createDimension("sigma", len(sigma))

# first argument is name of variable, 
# second is datatype,
# third is a tuple with the names of dimensions.

lat_netcdf = ncfile.createVariable("nav_lat", np.dtype("float32").char, ("y", "x"))
lat_netcdf[:] = nav_lat
lat_netcdf.units = "deg"
lat_netcdf.long_name = "latitude"

lon_netcdf = ncfile.createVariable("nav_lon", np.dtype("float32").char, ("y", "x"))
lon_netcdf[:] = nav_lon
lon_netcdf.units = "deg"
lon_netcdf.long_name = "longitude"

sig_netcdf = ncfile.createVariable(args.sigma_var, np.dtype("float32").char, ("sigma"))
sig_netcdf[:] = sigma
sig_netcdf.units = "kg/m3"
sig_netcdf.long_name = args.sigma_var

data1_netcdf = ncfile.createVariable("isodep_tave", np.dtype("float32").char, ("sigma", "y", "x"),
                                     fill_value = 1e20)
data1_netcdf[:] = isodep_tave
data1_netcdf.units = "m"

data2_netcdf = ncfile.createVariable("epe", np.dtype("float32").char, ("sigma", "y", "x"),
                                     fill_value = 1e20)
data2_netcdf[:] = epe
data2_netcdf.units = "m2/s2"

data3_netcdf = ncfile.createVariable("epe_zint", np.dtype("float32").char, ("y", "x"),
                                     fill_value = 1e20)
data3_netcdf[:] = epe_zint
data3_netcdf.units = "m3/s2"


print(" ")
print("returning final time-averaged field %s" % args.file_out)
