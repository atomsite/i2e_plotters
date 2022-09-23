#!/usr/bin/env python3

import simsetup
from numpy import logspace
import argparse
import sys


def get_abundance_range(abun):
  if abun == None:
    abun_range = [0.0]
  else:
    abun_min   = abun[0]
    abun_max   = abun[1]
    abun_step  = abun[2]
    abun_range = logspace(abun_min,abun_max,abun_step)
  return abun_range

def main(args):

  if args.log_range == True:
    print("Using log space abundance ranges!")
    if len(args.al26_abundances) != 3 or len(args.fe60_abundances) != 3:
      print("Abundances arguments must be formatted in the form <log10(min)> <log10(max)> <nsteps>")
      sys.exit()
    al26_range = get_abundance_range(args.al26_abundances)
    fe60_range = get_abundance_range(args.fe60_abundances)
  else:
    al26_range = args.al26_abundances
    fe60_range = args.fe60_abundances

  # Use a nested for loop to extract all permutations
  nsims = 0
  for planetoid_size in args.planetoid_sizes:
    for core_ratio in args.core_ratios:
      for al26_abun in al26_range:
        for fe60_abun in fe60_range:
          # Get a name for the working directory
          work_dir = "sim_p_{:.2f}_c_{:.2f}_al_{:.2f}_fe_{:.2f}".format(
            planetoid_size,
            core_ratio,
            al26_abun,
            fe60_abun)
          simsetup.gen_sim(planetoid_size,
                           core_ratio,
                           al26_abun,
                           fe60_abun,
                           args.exit_time,
                           args.init_time,
                           work_dir)
          simsetup.write_yaml(work_dir,planetoid_size,core_ratio,al26_abun,fe60_abun,args.exit_time,args.init_time)
          nsims += 1
  
  print("Configured {} simulations!".format(nsims))

if __name__ == "__main__":
  parse = argparse.ArgumentParser(description="Build a series of simulations to be processed using i2elvis")
  parse.add_argument("-s","--planetoid_sizes",nargs="+",help="A list containing the planetesimal sizes that are being used",type=float,default=[50.])
  parse.add_argument("-c","--core_ratios",nargs="+",help="A list containing the core ratios to use",type=float,default=[0.])
  parse.add_argument("-l","--log_range", default=False, action="store_true", help="Use -al and -fe flags to determine log spacing, in the form <log10(min)> <log10(max)> <steps>")
  parse.add_argument("-al","--al26_abundances",nargs="+",help="Al26 abundances relative to earth",type=float,default=None)
  parse.add_argument("-fe","--fe60_abundances",nargs="+",help="Fe60 abundances relative to earth",type=float,default=None)
  parse.add_argument("-i","--init_time",type=float,help="Initialisation time (yr)",default=1.0e6)
  parse.add_argument("-f","--exit_time",type=float,help="Initialisation time (yr)",default=11.0e6)
  args = parse.parse_args()
  main(args)