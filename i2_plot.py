from yaml import safe_load
from pandas import read_csv
import matplotlib.pyplot as plt

class Simulation:
  def __init__(self,path) -> None:
    self.path = path
    self.data_filename = path + "/hydrous_silicates.t3c"
    self.head_filename = path + "/simulation.yaml"
    # Import data from hydrous_silicates file
    self.read_data()
    # Import simulation info from simulation.yaml
    self.read_header()
    self.time = self.data.time - (self.sim["init_time"]/1e6)
  def read_data(self):
    self.data = read_csv(self.data_filename,sep=" ")
    self.data.columns = ["time",
                         "sol_frac",
                         "liq_frac",
                         "hydrous_frac",
                         "primitive_frac",
                         "n2co2_frac",
                         "cocl_frac",
                         "h2o_frac",
                         "phyllo1_frac",
                         "phyllo2_frac",
                         "phyllo3_frac",
                         "phyllo4_frac",
                         "perco_frac",
                         "melt1_frac",
                         "melt2_frac",
                         "maxtk",
                         "t_max_body",
                         "meantk",
                         "t_mean_body",
                         "count_toohot"]
  def read_header(self):
    with open(self.head_filename, "r") as file:
      self.sim = safe_load(file)
      self.header = self.sim # Whoops


def plot(fig,sim,parameter,label=None,col=None,linestyle="solid",logy=True,alpha=1.0):
  if col == None:
    if logy == True:
      plt.loglog(sim.time,sim.data[parameter],label=label,linestyle=linestyle,alpha=alpha)
    else:
      plt.semilogx(sim.time,sim.data[parameter],label=label,linestyle=linestyle,alpha=alpha)
  else:
    col = "C{}".format(col)
    if logy == True:
      plt.loglog(sim.time,sim.data[parameter],col,label=label,linestyle=linestyle,alpha=alpha)
    else:
      plt.semilogx(sim.time,sim.data[parameter],label=label,linestyle=linestyle,alpha=alpha)

def init_plot(width=5,height=5,grid=True):
  plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
  })
  # plt.style.use('seaborn-darkgrid')
  fig = plt.figure(figsize=(width,height))
  if grid == True:
    plt.grid(True, which="both", ls="dotted",alpha=0.3)
  return fig