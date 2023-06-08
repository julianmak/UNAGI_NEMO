#!/usr/bin/env python3

# JM: 26 Aug 2022
# script to generate the depth-integrated eke
# do this from the grid_U files directly

import numpy as np
import glob, netCDF4, copy, sys

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = 
    """
    Script to generate the depth-integrated eke.
    
    Computes <u>^2 and <u^2> from grid_U files, assumes the files are downsized
    and small enough that one could just brute force load everything (otherwise
    go in the code and load per layer say)
    """)
    
# fixed arguments
parser.add_argument("data_dir", type = str,
                    help = "specify data directory")
parser.add_argument("keys_U",     type = str, 
                    help = "specify string in filename to grab (e.g. \"*grid_U*\", note the quote marks, don't include the .nc extension)")
parser.add_argument("keys_V",     type = str, 
                    help = "specify string in filename to grab (e.g. \"*grid_V*\", note the quote marks, don't include the .nc extension)")

# optional arguments
parser.add_argument("--file_out", type = str, 
                    help = "specify output name (e.g. eke_zint_tave.nc)",
                    default = "eke_zint_tave.nc")
parser.add_argument("--lquery",
                    help = "print out the variables available",
                    action = "store_true")
parser.add_argument("--u_var", type = str,
                    help = "variable name for u (default = uoce)",
                    default = "uoce")
parser.add_argument("--v_var", type = str,
                    help = "variable name for v (default = voce)",
                    default = "voce")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------

# initialisation

# grab file names  
filenames = args.data_dir + args.keys_U + ".nc"
file_U_list = []
for file in glob.glob(filenames):
  file_U_list.append(file)
    
if len(file_U_list) == 0:
  sys.exit("no U files grabbed in %s, are you in the right folder?" % filenames)
  
filenames = args.data_dir + args.keys_V + ".nc"
file_V_list = []
for file in glob.glob(filenames):
  file_V_list.append(file)
  
if len(file_V_list) == 0:
  sys.exit("no V files grabbed in %s, are you in the right folder?" % filenames)
  
if len(file_V_list) != len(file_U_list):
  print(f"U_file list = {len(file_U_list)}")
  print(f"V_file list = {len(file_V_list)}")
  sys.exit("U and V files length mismatch, stopping..." % filenames)

# sort it according to the timestamps
file_U_list.sort()
file_V_list.sort()

#--------------------------------------------------------
# cycle through the file lists and compile the moc
print("%g files found, cycling through them..." % len(file_U_list))

for i in range(len(file_U_list)):

  # load e3t and masks here
  # if the filesize is too big, load it per layer say
  data = netCDF4.Dataset(args.data_dir + "mesh_mask.nc")
  e3t = data.variables["e3t_0"][0, :, :, :]
  e3u = data.variables["e3u_0"][0, :, :, :]
  e3v = data.variables["e3v_0"][0, :, :, :]
  e1u, e2u = data.variables["e1u"][0, :, :], data.variables["e2u"][0, :, :]
  e1v, e2v = data.variables["e1v"][0, :, :], data.variables["e2v"][0, :, :]
  e1t, e2t = data.variables["e1t"][0, :, :], data.variables["e2t"][0, :, :]
  tmask = data.variables["tmask"][0, :, :, :]
  
  nav_lon  = data.variables["nav_lon"][:, :]
  nav_lat  = data.variables["nav_lat"][:, :]
  
  try:
    npk = data.dimensions["z"].size
  except KeyError:
    print("default key \"z\" for npk not found, trying another possibility \"nav_lev\"")
    npk = data.dimensions["nav_lev"].size
  npj = data.dimensions["y"].size
  npi = data.dimensions["x"].size
  data.close()
  
  e1e2t = e1t * e2t
  e1e2u = e1u * e2u
  e1e2v = e1v * e2v
  
  del e1u, e2u, e1v, e2v, e1t, e2t
  
  print(" ")
  print("working in file = %g / %g" % (i + 1, len(file_U_list)))
  data_U = netCDF4.Dataset(file_U_list[i])
  data_V = netCDF4.Dataset(file_V_list[i])
  npt = data_U.dimensions["time_counter"].size
  if i == 0:
    u_tavg , v_tavg  = np.zeros((npk, npj, npi)), np.zeros((npk, npj, npi))
    u2_tavg, v2_tavg = np.zeros((npk, npj, npi)), np.zeros((npk, npj, npi))
    print("  found npk = %i, npj = %i, npi = %i, npt = %i" % (npk, npj, npi, npt))
    avg_total = len(file_U_list) * npt
    print("  going to average over npt * len(file_list) = %i records" % avg_total)

  if args.lquery:
    for name, variable in data.variables.items():
      for attrname in variable.ncattrs():
        if attrname == "standard_name":
          print("{} -- {}".format(name, getattr(variable, attrname)))

    data.close()

    print(" ")

    sys.exit("finished query, exiting gen_ekezint_tave...")
    
  else:
  
    # I/O intensive way of loading
    ##########################################################
    # !!! Consider using xarray for these kind of things...
    ##########################################################
      
    for kt in range(npt):
      if kt % 10 == 0:
        print(f"working in time kt = {kt+1} / {npt}")
      u_tavg  += data_U.variables[args.u_var][kt, :, :, :]    / avg_total
      u2_tavg += data_U.variables[args.u_var][kt, :, :, :]**2 / avg_total
      
      v_tavg  += data_V.variables[args.v_var][kt, :, :, :]    / avg_total
      v2_tavg += data_V.variables[args.v_var][kt, :, :, :]**2 / avg_total

  data_U.close()
  data_V.close()
  
#--------------------------------------------------------
# compute eke as <u^2> - <u>^2

print("computing eke...")

ue = u2_tavg - (u_tavg)**2
ve = v2_tavg - (v_tavg)**2

eke = np.zeros((npk, npj, npi))

for k in range(npk-1):
  zztmp = 0.25 / e1e2t / e3t[k, :, :]
  eke[k, 1:npj, 1:npi] = ( zztmp[1:npj, 1:npi] 
                         * (  ue[k, 1:npj  , 0:npi-1] * e1e2u[1:npj  , 0:npi-1] * e3u[k, 1:npj  , 0:npi-1]
                            + ue[k, 1:npj  , 1:npi  ] * e1e2u[1:npj  , 1:npi  ] * e3u[k, 1:npj  , 1:npi  ]
                            + ve[k, 0:npj-1, 1:npi  ] * e1e2v[0:npj-1, 1:npi  ] * e3v[k, 0:npj-1, 1:npi  ]
                            + ve[k, 1:npj  , 1:npi  ] * e1e2v[1:npj  , 1:npi  ] * e3v[k, 1:npj  , 1:npi  ]
                           )
                         )

ekezint  = np.sum(eke * tmask * e3t, axis=0)

#--------------------------------------------------------
# write the file

print("writing...")

# open a new netCDF file for writing.
ncfile = netCDF4.Dataset(args.data_dir + args.file_out, "w", format = "NETCDF4") 
ncfile.title = "(specific) eke depth-integrated file"

# create the dimensions.
ncfile.createDimension("x", nav_lon.shape[1])
ncfile.createDimension("y", nav_lon.shape[0])

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

data_netcdf = ncfile.createVariable("eke_zint", np.dtype("float32").char, ("y", "x"),
                                    fill_value = 1e20)
data_netcdf[:] = ekezint
data_netcdf.units = "m3/s^2"


print(" ")
print("returning final time-averaged field %s" % args.file_out)
