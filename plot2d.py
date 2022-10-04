# -*- coding: utf-8 -*-
"""
Visualisation routine for I2ELVIS_planet .prn output files
originally based on the Matlab 2D plotting routine by Gregor Golabek (GJG)
Created on 2015-05-15 by Tim Lichtenberg (TL)
Edit 2015-05-15: Change contour to heatmap, various layout changes
Edit 2018-04-20: TL: Added changes to adapt to github usage
Edit 2022-09-11: JE: Hasty functionalisation, run individual inputs for gnu parallel direct in folder
"""

# Import packages
import matplotlib # Plotting routine
matplotlib.use('Agg') # to counteract problems with torque usage
import os #https://stackoverflow.com/questions/3964681/find-all-files-in-directory-with-extension-txt-with-python
import sys
import glob
import struct # Decode binary data
import numpy as np 
from natsort import natsorted #https://pypi.python.org/pypi/natsort
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, DrawingArea, HPacker
from matplotlib import rc
from matplotlib import cm
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
# import colorlover as cl #http://nbviewer.ipython.org/gist/jackparmer/4696226c9b6b931bbdf6
#from IPython.display import HTML
import locale # http://stackoverflow.com/questions/1823058/how-to-print-number-with-commas-as-thousands-separators
import datetime
import argparse
# from subprocess import check_output
# from uuid import uuid4
import gzip

class Image:
    def __init__(self,filename) -> None:
        self.filename = filename
        self.read_data(self.filename)
        return
    class Var:
        def __init__(self,var_name,symbol,units,data) -> None:
            self.var_name = var_name
            self.symbol   = symbol
            self.units    = units
            self.data     = data
            self.label    = "{} ({})".format(symbol,units)
            # Calculate statistics
            self.mean     = np.mean(data)
            # Calculate quartiles
            self.quart = calc_array_quartiles(data)
            # Calculate other values, trivial calculation time
            self.min    = self.quart[0]
            self.max    = self.quart[4] 
            self.median = self.quart[2] # Median, synonym
            self.iqr    = self.quart[3] - self.quart[1] # Interquartile range
    def read_data(self,filename):
        """
        Subset of Tim's plotting code that just reads in the data from the simulation file
        """
        if filename[-3:] == ".gz":
            # print("Reading compressed data...")
            fdata = gzip.open(filename, "rb")
        else:
            # print("Reading data...")
            fdata = open(filename, "rb")

        # Read sizes of variables
        # A = struct.unpack("BBBB",fdata.read(4))[0]
        A           = fdata.read(4)
        # Read model parameters
        # Grid resolution
        xnumx       = struct.unpack("q",fdata.read(8))[0]
        ynumy       = struct.unpack("q",fdata.read(8))[0]
        # Markers per cell
        mnumx       = struct.unpack("q",fdata.read(8))[0]
        mnumy       = struct.unpack("q",fdata.read(8))[0]
        # Number of markers
        marknum     = struct.unpack("q",fdata.read(8))[0]
        # Model sizes
        xsize       = struct.unpack("d",fdata.read(8))[0]
        ysize       = struct.unpack("d",fdata.read(8))[0]
        # Efficieny of impact heating 
        gamma_eff   = struct.unpack("d",fdata.read(8))[0]
        # Thermal memory of iron material
        memory_fe   = struct.unpack("d",fdata.read(8))[0]
        # Thermal memory of silicate material
        memory_si   = struct.unpack("d",fdata.read(8))[0]
        # Porosity parameter
        por_init    = struct.unpack("d",fdata.read(8))[0]
        # Porosity parameter
        growth_model = struct.unpack("q",fdata.read(8))[0]
        # Grain size parameter
        gr_init     = struct.unpack("d",fdata.read(8))[0]
        # Grid resolution in z direction
        znumz       = struct.unpack("q",fdata.read(8))[0]
        # 2D cylindrical or 3D spherical
        corr2d3d    = struct.unpack("d",fdata.read(8))[0]
        # Pressure value
        pinit       = fdata.read(40)
        # Gravity
        GXKOEF      = struct.unpack("d",fdata.read(8))[0]
        GYKOEF      = struct.unpack("d",fdata.read(8))[0]
        # Ambient temperature for sticky air
        tmp_ambient = struct.unpack("d",fdata.read(8))[0]
        # Time to exit simulation
        timeexit    = struct.unpack("d",fdata.read(8))[0]
        # 26Al/27Al ratio
        al2627_init = struct.unpack("d",fdata.read(8))[0]
        # 60Fe/56Fe ratio
        fe6056_init = struct.unpack("d",fdata.read(8))[0]
        # Number of rocks
        rocknum     = struct.unpack("i",fdata.read(4))[0]
        # Number of Boundary conditions
        bondnum     = struct.unpack("q",fdata.read(8))[0]
        # Stage
        n1          = struct.unpack("i",fdata.read(4))[0]
        # Time
        Time        = struct.unpack("d",fdata.read(8))[0]
        
        # print A, xnumx, ynumy, mnumx, mnumy, marknum
        # print xsize, ysize, gamma_eff, memory_fe, memory_si
        # print por_init, znumz, corr2d3d, pinit, GXKOEF
        # print GYKOEF, rocknum, bondnum, n1, Time

        # Define cursor positions within the binary file in bytes
        # 1*int (A),              x*int,    y*float/long
        curpos_nodes    = 4     + 2*4   +   28*8 + rocknum*(8*24+4)
        # Skip all nodes
        curpos_grid     = curpos_nodes+(4*20+8*4)*xnumx*ynumy
        # Skip gridlines and boundary conditions
        curpos_markers  = curpos_grid+4*(xnumx+ynumy)+(4*4+8*3)*(bondnum-1)

        # print ("Read in node properties...")

        # Skip rock properties
        #curpos0=4+2*4+22*8+rocknum*(8*24+4)
        fdata.seek(curpos_nodes)

        prs   = np.zeros((ynumy,xnumx))  # pressure
        vx    = np.zeros((ynumy,xnumx))  # x component of the velocity
        vy    = np.zeros((ynumy,xnumx))  # y component of the velocity
        exx   = np.zeros((ynumy,xnumx))  # xx component of strain rate tensor
        eyy   = np.zeros((ynumy,xnumx))  # yy component of strain rate tensor
        exy   = np.zeros((ynumy,xnumx))  # xy component of strain rate tensor
        sxx   = np.zeros((ynumy,xnumx))  # xx component of stress tensor
        syy   = np.zeros((ynumy,xnumx))  # yy component of stress tensor
        sxy   = np.zeros((ynumy,xnumx))  # xy component of stress tensor
        rho   = np.zeros((ynumy,xnumx))  # density
        eta   = np.zeros((ynumy,xnumx))
        etd   = np.zeros((ynumy,xnumx))
        mu    = np.zeros((ynumy,xnumx))
        gp    = np.zeros((ynumy,xnumx))  # gravity potential
        ep    = np.zeros((ynumy,xnumx))  # viscosity 1
        et    = np.zeros((ynumy,xnumx))  # viscosity 2
        tmp   = np.zeros((ynumy,xnumx))  # temperature
        cp    = np.zeros((ynumy,xnumx))  # heat capacity
        kt    = np.zeros((ynumy,xnumx))  # thermal conductivity
        ht    = np.zeros((ynumy,xnumx))  # radiogenic heating
        eii   = np.zeros((ynumy,xnumx))  # strain rate tensor
        sii   = np.zeros((ynumy,xnumx))  # stress tensor
        dis   = np.zeros((ynumy,xnumx))  # dissipation
        fm    = np.zeros((ynumy,xnumx))  # silicate melt fraction
        fm_fe = np.zeros((ynumy,xnumx))  # iron melt fraction
        radius= np.zeros((ynumy,xnumx))  # distance from planetary center [m]

        # Read nodes information
        for i in range(0,xnumx):
            for j in range(0,ynumy):
                prs[j,i]=struct.unpack("f",fdata.read(4))[0]
                vx[j,i]=struct.unpack("f",fdata.read(4))[0]
                vy[j,i]=struct.unpack("f",fdata.read(4))[0]
                vbuf1=fdata.read(3*8)
                exx[j,i]=struct.unpack("f",fdata.read(4))[0]
                eyy[j,i]=struct.unpack("f",fdata.read(4))[0]
                exy[j,i]=struct.unpack("f",fdata.read(4))[0]
                sxx[j,i]=struct.unpack("f",fdata.read(4))[0]
                syy[j,i]=struct.unpack("f",fdata.read(4))[0]
                sxy[j,i]=struct.unpack("f",fdata.read(4))[0]
                rho[j,i]=struct.unpack("f",fdata.read(4))[0]
                eta[j,i]=struct.unpack("f",fdata.read(4))[0]
                etd[j,i]=struct.unpack("f",fdata.read(4))[0]
                mu[j,i]=struct.unpack("f",fdata.read(4))[0]
                gp[j,i]=struct.unpack("f",fdata.read(4))[0]
                ep[j,i]=struct.unpack("f",fdata.read(4))[0]
                et[j,i]=struct.unpack("f",fdata.read(4))[0]
                tmp[j,i]=struct.unpack("f",fdata.read(4))[0]
                vbuf3=fdata.read(1*8)
                cp[j,i]=struct.unpack("f",fdata.read(4))[0]
                kt[j,i]=struct.unpack("f",fdata.read(4))[0]
                ht[j,i]=struct.unpack("f",fdata.read(4))[0]
        # Compute distance of each grid point from planetary center [m]
        for i in range(1,xnumx):
            for j in range(1,ynumy):    
                radius[j,i]=(((i-1)*(xsize/(xnumx-1))-(xsize/2.000))**(2.000)+((j-1)*(ysize/(ynumy-1))-(ysize/2.000))**(2.000))**(0.500)
        # Joe edit: Extract variables
        self.xnumx, self.ynumy = xnumx, ynumy
        self.xsize, self.ysize = xsize, ysize
        # Create extent for plotting centred at midpoint
        xx = xsize / 2000 # divide by 2, convert to km
        yy = ysize / 2000 # divide by 2, convert to km
        xrange = [-xx,xx]
        yrange = [-yy,yy]
        self.extent = [xrange[0],xrange[1],yrange[0],yrange[1]]

        # Get spacing of each cell
        self.x_space = xsize / xnumx # X direction
        self.y_space = ysize / ynumy # Y direction, This should be the same as X dir
        # Array starting positions
        self.x_start = -xx + (self.x_space / 2.0)
        self.y_start = -yy + (self.y_space / 2.0)
        # Get other parameters
        # Abundance parameters
        al_canonical = 5.25e-5 # Canonical solar system abundance of 26Al
        fe_canonical = 1.15e-8 # Canonical solar system abundance of 60Fe
        self.alratio = al2627_init
        self.feratio = fe6056_init
        self.al_solar_ratio = al2627_init / al_canonical
        self.fe_solar_ratio = fe6056_init / fe_canonical
        # Time parameters
        self.t       = Time
        self.t_myr   = Time / 1e6
        self.t_exit  = timeexit
        self.n_out   = n1
        
        # Joe edit: Extract data from cells
        self.rho=self.Var("Density",
                          r'$\rho',
                          r'$\mathrm{kg} \, \mathrm{m}^{-3}$',
                          rho)
        self.temp = self.Var("Temperature",
                              r"T",
                              r"K",
                              tmp)
        self.kt = self.Var("Thermal conductivity",
                           r'$k_{\mathrm{eff}}',
                           r'\mathrm{W} \, \mathrm{m}^{-1} \, \mathrm{K}^{-1}$',
                           kt)
        self.press = self.Var("Pressure",
                              r"p",
                              r"Pa",
                              prs)
        self.radius = self.Var("Radius",
                               r"$\mathrm{r}$",
                               r"$\mathrm{m}",
                               radius)
        self.ht = self.Var("Radiogenic heating",
                           r"$\mathcal{Q}$",
                           r"$\mathrm{w}\,\mathrm{kg}^{-1}",
                           ht)
    
    def planetesimal_radii(self,mantle_radius,core_mantle_ratio):
        self.mantle_radius = mantle_radius
        self.core_mantle_ratio = core_mantle_ratio
        self.core_radius = mantle_radius * core_mantle_ratio

    def get_temp_statistics(self):
        try:
            self.mantle_radius
            self.core_radius
        except:
            raise("Missing radius parameters of simulation, run planetesimal_radii(mantle_radius,core_mantle_ratio)")
        try: 
            self.core_temp
            self.mantle_temp
        except:
            self.calc_core_temps()
            self.calc_mantle_temps()
        self.core_mean_temp    = np.mean(self.core_temp)
        self.core_min_temp     = np.min(self.core_temp)
        self.core_max_temp     = np.max(self.core_temp)
        self.core_stdev_temp   = np.std(self.core_temp)
        self.mantle_mean_temp  = np.mean(self.mantle_temp)
        self.mantle_min_temp   = np.min(self.mantle_temp)
        self.mantle_max_temp   = np.max(self.mantle_temp)
        self.mantle_stdev_temp = np.std(self.mantle_temp)

    def calc_region_temps(self):
        self.calc_core_temps()
        self.calc_mantle_temps()

    def calc_core_temps(self):
        try:
            test = self.mantle_radius
            test = self.core_radius
        except:
            raise("Missing radius parameters of simulation, run planetesimal_radii(mantle_radius,core_mantle_ratio)")
        self.core_temp = radius_subset(self.temp.data,self.radius.data,0.0,self.core_radius,padding=0.05)
        # self.core_temp = gut_ambient(self.core_temp)
    
    def calc_mantle_temps(self):
        try:
            test = self.mantle_radius
            test = self.core_radius
        except:
            raise("Missing radius parameters of simulation, run planetesimal_radii(mantle_radius,core_mantle_ratio)")
        self.mantle_temp = radius_subset(self.temp.data,self.radius.data,self.core_radius,self.mantle_radius,padding=0.05)
        # self.mantle_temp = gut_ambient(self.mantle_temp)

    def calc_desiccation(self):
        # np.zeros((self.ynumy,self.xnumx),dtype=np.int8)
        # for i in range(self.ynumy):
        #     for j in range(self.xnumx):
        #         if self.temp.data[i,j] >= 



        return

    def calc_temperature_ranges(self):
        

        # Get temperature ranges from 

        return

def calc_array_quartiles(data):
    quart = [] # Make empty list
    quart.append(np.min(data))
    quart.append(np.quantile(data,0.25)) # 1st quartile
    quart.append(np.quantile(data,0.50)) # 2nd quartile
    quart.append(np.quantile(data,0.75)) # 3rd quartile
    quart.append(np.max(data))
    return quart

def radius_subset(quant,radius,r_min,r_max,padding=0.0):
    q = []
    del_r = r_max - r_min
    r_pad = del_r * padding
    r_min = r_min + r_pad
    r_max = r_max - r_pad
    ny,nx = np.shape(quant)
    for i in range(ny):
        for j in range(nx):
            rad = radius[i,j]
            if rad > r_min and rad < r_max:
                q.append(quant[i,j])
    return np.asarray(q)

def gut_ambient(temps):
    t_ambient = 155.
    new_temps = []
    for temp in temps:
        if temp > t_ambient:
            new_temps.append(temp)
    return np.asarray(new_temps)

def plot_circles(ax,sim,colour="r",alpha=0.3):
    """
    Plot circles
    """
    def draw(r):
        theta = np.linspace(0,2*np.pi,1000)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        ax.plot(x,y,color=colour,alpha=alpha)
    draw(sim.header["core_radius"])
    draw(sim.header["mantle_radius"])
    return ax

def plot_split(ax,data,col):
    """
    Plot split lines of temporal data
    Plots lines for:
        - min
        - max
        - mean
        - +- stdev
    Plots fills for:
        - min->max
        - -stdev->+stdev
    Inputs:
        - ax: input figure
        - data pandas array generated from csv from postprocess.py
        - col: plot colour, either numeric ("CN") or specific ("r","k",etc.)
    Outputs:
        - ax: parts of a plot
    """
    if type(col) == int:
        col="C"+str(col)
    ax.fill_between(data.t,data["min"],data["max"],color=col,alpha=0.25)
    ax.fill_between(data.t,data["mean"]-data["std"],data["mean"]+data["std"],color=col,alpha=0.25)
    ax.plot(data.t,data["min"],color=col)
    ax.plot(data.t,data["max"],color=col)
    ax.plot(data.t,data["mean"],color=col)
    ax.plot(data.t,data["mean"]-data["std"],color=col)
    ax.plot(data.t,data["mean"]+data["std"],color=col)
    return ax

class Desiccation:
    def __init__(self,folder_name) -> None:
        """
        Calculate desiccation levels
        Also includes plot routines
        """
        filename = folder_name+"/desiccation_times.npz"
        try:
            data = np.load(filename)
            self.vap  = data["vap"]
            self.melt = data["melt"]
        except:
            print("Something's up with reading!")
            raise
    def generate_img(self,t):
        """
        For a given time, generate a numpy array
        Key:
            - If cell is frozen    = 0
            - If cell is melting   = 1
            - If cell is vaporised = 2
        """
        res = np.shape(self.vap)
        img = np.zeros(res)
        for i in range(res[0]):
            for j in range(res[1]):
                vap_time  = self.vap[i,j]
                melt_time = self.melt[i,j]
                pix = 0
                if t >= melt_time:
                    pix = 1
                if t >= vap_time:
                    pix = 2
                img[i,j] = pix
        # Store most recent image array in class
        self.last_img = img
        # Also return if needed
        return img
    def plot(self,ax):
        ax.imshow(self.last_img)
        return ax
    def desiccation_percent(self):
        return
    def desiccation_raw(self):
        """
        Get the actual number of dessicated cells
        """
        return
