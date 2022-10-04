
import numpy as np
from plot2d import Image,Desiccation,plot_split
from i2_plot import Simulation
from glob import glob
from tqdm import tqdm
import pandas as pd

def desiccation_times(folder):
  # Constants
  t_vap  = 1223.15 # Full boiling temperature (K)
  t_melt = 273.15  # Melting temperature (K)
  # Read simulation parameters
  sim = Simulation(folder)
  mantle_radius = sim.header["mantle_radius"] * 1000
  core_ratio = sim.header["core_fraction"]
  init_time = sim.header["init_time"]
  # Determine how many snapshots there are
  files = sorted(glob(folder+"/*.prn*"))
  if len(files) == 0:
    data_folder = folder+"/prnfiles/"
    files  = sorted(glob(data_folder+"*.prn.gz"))
    nfiles = len(files)
    img0   = Image(data_folder+"dd_0.prn.gz")
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
  # Create pandas arrays to handle temperature data
  # Pandas is being used for columns info, csv writing
  # Numpy binary blobs used for mesh data as it is faster and more useful for use-case

  filename = data_folder+"/dd_0.prn.gz"
  image = Image(filename)
  image.planetesimal_radii(mantle_radius,core_ratio)
  image.get_temp_statistics()

  t = image.t - init_time
  df_mant = describe(t,image.mantle_temp)
  df_core = describe(t,image.core_temp)

  # Get values, assume that no heating occurs at t=0
  for n in tqdm(range(1,nfiles)):
    filename = data_folder+"dd_"+str(n)+".prn.gz"
    image = Image(filename)
    t = image.t - init_time
    # Calculate core and mantle temperatures, a very small degree of wiggle room is used
    image.planetesimal_radii(mantle_radius,core_ratio) # Code needs mantle and core radii
    image.calc_region_temps() # 
    # Write mantle values
    df_mant = pd.concat([df_mant,describe(t,image.mantle_temp)],join="inner",ignore_index=True)
    df_core = pd.concat([df_core,describe(t,image.core_temp)],join="inner",ignore_index=True)
    # Write desiccation data
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
  np.savez_compressed(folder+"/desiccation_times",vap=vap_arr,melt=melt_arr)
  # Write temperature statistics to csv files
  df_mant.to_csv(folder+"/mantle_temps.csv")
  df_core.to_csv(folder+"/core_temps.csv")
  # Finished! Returns arrays if used in a function
  return vap_arr, melt_arr

def plot_test(folder):
  df_core = pd.read_csv(folder+"/core_temps.csv")
  df_mant = pd.read_csv(folder+"/mantle_temps.csv")
  img = Desiccation(folder)
  img.generate_img(10e6)

  import matplotlib.pyplot as plt
  f, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 1, 1]})
  a0 = img.plot(a0)
  a1 = plot_split(a1,df_core,0)
  a2 = plot_split(a2,df_mant,1)
  a1.set_xscale("log")
  a2.set_xscale("log")
  plt.savefig("test.pdf")
  
def describe(t,array):
  df = pd.DataFrame(array)
  dd = df.describe().T
  dd["t"] = t
  return dd

if __name__ == "__main__":
  import sys
  folder = sys.argv[1]
  # desiccation_times(folder)
  plot_test(folder)