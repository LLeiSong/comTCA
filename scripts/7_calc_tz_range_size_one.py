# Title     : Calculate species-level richness index
# Objective : The main function to calculate the species-level richness index
# Created by: Lei Song
# Created on: 07/26/23
# The downstream for calc_tz_range_size_scheduler.R

# Load libraries
import os
import re
import glob
import click
from os.path import exists, join, basename, isfile
from osgeo import gdal
import numpy as np
import multiprocessing as mp
import pandas as pd
from sys import exit

# Define the function to calculate the range size
def get_ranges_size(fname):
  img = gdal.Open(fname)
  band = img.GetRasterBand(1)
  values = band.ReadAsArray()
  
  return values.sum()

@click.command()
@click.option('--fname', type=str,
              help="The path of the file to process.")
def main(fname):
  print(fname)
  
  # Set directories
  ref_rg_dir = join("/scratch/lsong36/comTCA/data/expert_range_maps",
                    "refined_range_rasters")
  tz_dir = join(ref_rg_dir, "tanzania")
  dst_dir = '/scratch/lsong36/comTCA/data/biodiveristy'
  if not exists(dst_dir):
    os.mkdir(dst_dir)
    
  # Define file to take results
  csv_path = join(dst_dir, "tanzania_range_size.csv")
  if not isfile(csv_path):
    exit("No csv file found!")
  
  # Calculate
  print("Calculate the value.")
  range_size = get_ranges_size(fname)
  species = basename(fname).replace(".tif", "").replace("_", " ")
  range_tbl = pd.DataFrame({"species": [species], "tanzania_range_size": [range_size]})
  print("Save out the result.")
  range_tbl.to_csv(csv_path, mode='a', index=False, header=False)
  
if __name__ == '__main__':
    main()
