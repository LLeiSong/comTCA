# Load libraries
import os
import re
import glob
from os.path import exists, join, basename
from osgeo import gdal
import numpy as np
import multiprocessing as mp
import pandas as pd

# Set directories
ref_rg_dir = join("/scratch/lsong36/comTCA/data/expert_range_maps",
                  "refined_range_rasters")
outside_dir = join(ref_rg_dir, "outside")
dst_dir = '/scratch/lsong36/comTCA/data/biodiveristy'
if not exists(dst_dir):
  os.mkdir(dst_dir)

# Define the function to calculate the range size
def get_ranges_size(fname):
  img = gdal.Open(fname)
  band = img.GetRasterBand(1)
  values = band.ReadAsArray()
  
  return values.sum()

# Get range sizes for the global refined maps
taxon_group = ["mammals", "amphibians", "birds", "reptiles"]
fnames = []
for taxon in taxon_group:
  fns = glob.glob(join(outside_dir, taxon, "*.tif"), recursive=True)
  if taxon == "birds":
    fns = list(filter(lambda fname: not re.search("_[NBR]{1}", fname), fns))
  fnames.extend(fns)

# Define file to take results
csv_path = join(dst_dir, "outside_range_size.csv")

# Parallel
species = [basename(fname).replace(".tif", "").replace("_", " ") for fname in fnames]
nthreads = 10
with mp.Pool(processes=nthreads) as pool:
  ranges = pool.map(get_ranges_size, fnames)
range_tbl = pd.DataFrame({"species": species,
                          "outside_range_size": ranges})
range_tbl.to_csv(csv_path, index=False)
