import sys
# sys.path.append("/local/jweatson/i2e_plotters")
from plot2d import Image, plot_circles
from i2_plot import Simulation
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import os
plt.rcParams.update({
  "text.usetex": True,
  "font.family": "serif",
  "font.serif": ["Computer Modern Roman"],
})

working_dir = "./"
filename = working_dir+sys.argv[1]
sim = Simulation(working_dir)
image = Image(filename)

norm = Normalize(vmin=0, vmax=1600)

f, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [3, 1, 1]})

timestr = "{:.3f}".format(image.t_myr)
plstr   = "{:.0f}".format(sim.header["mantle_radius"])
phistr  = "{:.2f}".format(sim.header["core_fraction"])
ironstr = "{:.0f}".format(sim.header["fe_ss_ratio"])

f.set_figheight(10)
f.set_figwidth(6)
pcm = a0.imshow(image.temp.data,extent=image.extent,norm=norm,cmap="inferno")
# Plot circles for mantle and core
a0.set_title(r"r$_p$ = {}$\,$km, $\Phi$ = {}, Fe = {}$\odot$, t = {} Myr".format(plstr,phistr,ironstr,timestr))
plot_circles(a0,sim)
a0.set_xlabel("X (km)")
a0.set_ylabel("Y (km)")
# Build the colour bar
cbar = f.colorbar(pcm,ax=a0,label=image.temp.label)

a1.plot(sim.time,sim.data["h2o_frac"],label="Water")
a1.plot(sim.time,sim.data["primitive_frac"],label="Primitive")
a1.plot(sim.time,sim.data["hydrous_frac"],label="Hydrous")
a1.axvline(image.t_myr - (sim.header["init_time"] / 1e6),color="r",linestyle="--")
a1.set_ylim(0,1)
a1.set_xscale("log")
a1.set_ylabel("$f_t/f_0$")
a1.legend(loc=1)
a1.set_xlim(1e-1,1e2)

a2.plot(sim.time,sim.data["meantk"],label="T$_\mathrm{mean}$")
a2.plot(sim.time,sim.data["maxtk"],label="T$_\mathrm{max}$")
a2.axvline(image.t_myr - (sim.header["init_time"] / 1e6),color="r",linestyle="--")
a2.set_xscale("log")
a2.set_ylabel("T (K)")
a2.legend(loc=1)
a2.set_xlim(1e-1,1e2)
a2.set_xlabel("$t$ (Myr)")

if not os.path.exists("figures"):
    os.makedirs("figures")

# figname = "figures/fig_{}.png".format(str(image.n_out).zfill(4))
import re
nn = re.findall('\d+', filename)
fign=str(nn[-1]).zfill(4)
# print(fign)
# import sys
# sys.exit()
figname = "figures/fig_"+fign+".png"
plt.savefig(figname,dpi=600)
print("Finished "+filename+" -> "+figname)