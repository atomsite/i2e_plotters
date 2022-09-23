
import numpy as np
from plot2d import Image,radius_subset
from glob import glob
from tqdm import tqdm

def dessication_times(folder):
  # Constants
  t_vap  = 1223.15 # Full boiling temperature (K)
  t_melt = 273.15  # Melting temperature (K)
  # Determine how many snapshots there are
  files = sorted(glob(folder+"/*.prn*"))
  if len(files) == 0:
    folder = folder+"/prnfiles/"
    files  = sorted(glob(folder+"*.prn.gz"))
    nfiles = len(files)
    img0   = Image(folder+"dd_0.prn.gz")
  else:
    print("Needs to be a complete dataset from a finished simulation!")
    sys.exit()
  # If still empty somethings gone wrong
  if len(files) == 0:
    print("No files found!")
    sys.exit()
  # Get data from first array
  nx = img0.xnumx
  ny = img0.ynumy
  tf = img0.t_exit
  # Build an array
  vap_arr  = np.full((nx,ny),tf)
  melt_arr = np.full((nx,ny),tf)
  # Make a file to store temperature data
  

  # Get values, assume that no heating occurs at t=0
  for n in tqdm(range(1,nfiles)):
    filename = folder+"dd_"+str(n)+".prn.gz"
    image = Image(filename)
    t = image.t
    for i in range(ny):
      for j in range(nx):
        temp = image.temp.data[i,j]
        if temp >= t_melt:
          if melt_arr[i,j] == tf:
            melt_arr[i,j] = t
        if temp >= t_vap:
          if vap_arr[i,j] == tf:
            vap_arr[i,j] = t
  # Write out to compressed blob
  np.savez_compressed(folder+"/dessication_times",vap=vap_arr,melt=melt_arr)
  # Finished! Returns arrays if used in a function
  return vap_arr, melt_arr

if __name__ == "__main__":
  import sys
  folder = sys.argv[1]
  dessication_times(folder)