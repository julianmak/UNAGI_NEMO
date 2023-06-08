#!/usr/bin/env python3

# JM: 31 Oct 2019
# script to generate the relevant data for time mean isopycnal depths
# for temperature co-ordinates, with time entries in the file

import numpy as np
import glob, netCDF4, copy, sys
from pyCDFTOOLS.eos import sigmai_dep

from numba import jit, int32 # need explicit type for use within jit

#--------------------------------------------------------
# define the argument parser
import argparse

parser = argparse.ArgumentParser(description = 
    """
    Script to generate the relevant data for time mean isopycnal depths. 
    
    Requires already processed files of the isopycnal depths, assumes that:
    * temperature is the co-ordinate [easy hack to make generic but can't be bothered now...]
    * multiple entries, load per-entry
    * number of records * (dep ** 2) is not so huge that it is going to overflow in single precision
    """)
    
# fixed arguments
parser.add_argument("data_dir", type = str,
                    help = "specify data directory")
parser.add_argument("keys",     type = str, 
                    help = "specify string in filename to grab (e.g. \"*isopyc_dep*\", note the quote marks, don't include the .nc extension)")

# optional arguments
parser.add_argument("--file_out", type = str, 
                    help = "specify output name (e.g. isodep_tave.nc)",
                    default = "isodep_tave.nc")
parser.add_argument("--lquery",
                    help = "print out the variables available",
                    action = "store_true")
parser.add_argument("--t_var", type = str,
                    help = "variable name for t (default = toce)",
                    default = "toce")
#parser.add_argument("--s_var", type = str,
#                    help = "variable name for s (default = vosaline)",
#                    default = "vosaline")

# collect arguments
args = parser.parse_args()

#--------------------------------------------------------
# define some subfunctions

@jit(nopython = True)
def isodep_loop(depi, wdep, jk, gdep, tmask, ibin, npiglo, npjglo):
  
    # Do the binning and convert k into sigma
    for ji in range(1, npiglo - 1):
      for jj in range(1, npjglo - 1):
        # cycle through the indices and bin according to density classes
        ib = ibin[jj, ji] - 1 # python indexing
        depi[ib, jj, ji] += gdep[jk] * tmask[jk, jj, ji]
        wdep[ib, jj, ji] +=            tmask[jk, jj, ji]

    return (depi, wdep)
    
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
  
  # load the mask [could probably save some I/O on this but putting it somewhere else...]
  data = netCDF4.Dataset(args.data_dir + "mesh_mask.nc")
  tmask   = data.variables["tmask"][0, :, :, :]
  gdep    =-data.variables["gdept_1d"][0, :]
  data.close()
  
  # cycle through file entries
  data = netCDF4.Dataset(file_list[i])
  
  npkglo = data.dimensions["deptht"].size
  npjglo = data.dimensions["y"].size
  npiglo = data.dimensions["x"].size
  npt    = data.dimensions["time_counter"].size
  
  print("  found npkglo = %i, npjglo = %i, npiglo = %i, npt = %i" % (npkglo, npjglo, npiglo, npt))
  print("  averaging is going to be over counts of non-zero entries in isodep")
  
  if i == 0:
    nav_lon = data.variables["nav_lon"][:, :]
    nav_lat = data.variables["nav_lat"][:, :]
    # call it sigma but it's actually going to be temperature
    nbins, sigmin, sigstp = 300, 34.5, 0.01
    sigma = sigmin + (np.linspace(1, nbins, nbins) - 0.5) * sigstp
    print("min and max density bin = %.4f, %.4f" % (sigma[0], sigma[nbins - 1]) )
    print("min and max db          = %.4f, %.4f" % (np.amin(np.diff(sigma)), np.amax(np.diff(sigma))))
    
    wdep_tave     = np.zeros((nbins, npjglo, npiglo))
    isodep_tave   = np.zeros((nbins, npjglo, npiglo))
    isodepsq_tave = np.zeros((nbins, npjglo, npiglo))

  if args.lquery:
    for name, variable in data.variables.items():
      for attrname in variable.ncattrs():
        if attrname == "standard_name":
          print("{} -- {}".format(name, getattr(variable, attrname)))

    data.close()

    print(" ")

    sys.exit("finished query, exiting gen_isodeps_tave...")
    
  else:
  
    # load per time-slice and do averaging
    # assume here that number of records * (dep ** 2) is not so huge that it is
    #   going to overflow in single precision
    for jt in range(0, npt): 
      print("working on jt = %i/%i..." % (jt + 1, npt))
      PT      = data.variables[args.t_var][jt, :, :, :]
      SP      = 35 * np.ones(PT.shape)

      depi = np.zeros((nbins, npjglo, npiglo))
      wdep = np.zeros((nbins, npjglo, npiglo))

      print("churning through the rebinning loops")
      print(" ")
      print("skipping last layer because it is just zeros")

      for jk in range(0, npkglo - 1):
        print("  working on jk = %i/%i..." % (jk + 1, npkglo - 1))

        # generate the neutral density 
#        dens = PT[jk, :, :]
        dens = sigmai_dep(PT[jk, :, :], SP[jk, :, :], 2000)

        # the mask takes care of the densities at the bottom
        zttmp  = dens * tmask[jk, :, :]

        # find bin numbers
        # a) faster procedure but only for uniform sigma grids
        ibin = (zttmp - sigmin) / sigstp
        ibin = np.minimum( np.maximum(ibin.astype(int), 1), nbins)

        # b) slower procedure but for general grids, would really like to get rid of the loop...
    #     ibin = np.zeros(zttmp.shape).astype(int)
    #     for ji in range(1, npiglo - 1):
    #         for jj in range(1, npjglo - 1):
    #             if tmask[0, jj, ji]: # avoid land points to save a bit of time
    #                 ibin[jj, ji] = (np.abs(sigma - zttmp[jj, ji])).argmin()
    #     ibin = np.minimum(np.maximum(ibin, 1), nbins - 1)
      
        # Do the binning and convert k into sigma
        depi, wdep = isodep_loop(depi, wdep, jk, gdep, tmask, ibin, npiglo, npjglo)

        # generate the raw dep with gaps in
        depi = np.where(wdep > 0.0, depi / wdep, 0.0)

        # interpolate to fill out some interior values
#        isopyc_dep = depi
        isopyc_dep = np.zeros((nbins, npjglo, npiglo))

        for ji in range(1, npiglo - 1):
          for jj in range(1, npjglo - 1):
            hits = np.where(depi[:, jj, ji] < 0.0)[0]
    
            if len(hits) > 0: # only do interpolation if it is non-empty
                isopyc_dep[hits[0]:hits[-1]+1, jj, ji] = (
                      np.interp(sigma[hits[0]:hits[-1]+1], sigma[hits], depi[hits, jj, ji])
                      )
    
    
        # only generate weights where there is a non-zero record
        wdep_tave     += np.where(isopyc_dep != 0.0, 1, 0)
        isodep_tave   += isopyc_dep
        isodepsq_tave += isopyc_dep * isopyc_dep
      
  data.close()

# average by the number of weights
isodep_tave   = np.where(wdep_tave > 0, isodep_tave   / wdep_tave, 0.0)
isodepsq_tave = np.where(wdep_tave > 0, isodepsq_tave / wdep_tave, 0.0)

#--------------------------------------------------------
# write the file

print("writing...")

# open a new netCDF file for writing.
ncfile = netCDF4.Dataset(args.data_dir + args.file_out, "w", format = "NETCDF4") 
ncfile.title = "tave files for isopycnal depths and its square in density co-ordinates, from UNAGI_R010 data"

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

sig_netcdf = ncfile.createVariable("sigma", np.dtype("float32").char, ("sigma"))
sig_netcdf[:] = sigma
sig_netcdf.units = "kg m-3"
sig_netcdf.long_name = "sigma"

data1_netcdf = ncfile.createVariable("isodep_tave", np.dtype("float32").char, ("sigma", "y", "x"),
                                     fill_value = 1e20)
data1_netcdf[:] = isodep_tave
data1_netcdf.units = "m"

data2_netcdf = ncfile.createVariable("isodepsq_tave", np.dtype("float32").char, ("sigma", "y", "x"),
                                     fill_value = 1e20)
data2_netcdf[:] = isodepsq_tave
data2_netcdf.units = "m2"


print(" ")
print("returning final time-averaged field %s" % args.file_out)
