# Title     : Calculate taxon-level species richness
# Objective : The main function to calculate the taxon-level species richness
# Created by: Lei Song
# Created on: 09/09/23
# Usage: python scripts/calc_species_richness.py --txn birds

# Load libraries
import os
import re
import glob
import click
from os.path import exists, join, basename
from osgeo import gdal
import numpy as np
import pandas as pd
import multiprocessing as mp

# Define the function to calculate the range size
def species_richness(fname, species_taxon):
  # Get numbers
  sps_nm = basename(fname).replace(".tif", "").replace("_", " ")
  sps_info = species_taxon.query("species == '{}'".format(sps_nm))
  t_weight = sps_info.iloc[0]['threat_weight'] * sps_info.iloc[0]['local_weight']
  range_size = sps_info.iloc[0]['refined_range_size']
  
  # Read layer
  img = gdal.Open(fname)
  band = img.GetRasterBand(1)
  sr = band.ReadAsArray()
  wsr = sr * t_weight
  rwsr = sr / range_size
  wrwsr = rwsr * t_weight
  img = None
  band = None
  
  return [sr, wsr, rwsr, wrwsr]

@click.command()
@click.option('--txn', type=str,
              help="The taxon of the species.")

def main(txn):
  # Set directories
  data_dir = "/scratch/lsong36/comTCA/data"
  range_dir = join(data_dir, "expert_range_maps/refined_range_rasters/tanzania_100")
  bio_dir = "/scratch/lsong36/comTCA/data/biodiveristy"
  
  # Summarize a bit
  species_info = pd.read_csv(join(bio_dir, "species_info.csv"))
  species_info = species_info.query("tanzania_range_size > 0")
  species_info = species_info.query("refined_range_size > 0")
  local_species = species_info.query("outside_range_size == 0")
  
  # Set weights
  lc_species_taxon = local_species.query("taxon == '{}'".format(txn))['species'].tolist()
  species_taxon = species_info.query("taxon == '{}'".format(txn)).copy()
  species_taxon.loc[:, 'local_weight'] = np.where(species_taxon['species'].isin(lc_species_taxon), 1, 1/2)
  # Whatever, use an embarassing way to set weights
  species_taxon.loc[:, 'threat_weight'] = 0
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'CR', 1, species_taxon.loc[:, 'threat_weight'])
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'EN', 1/2, species_taxon.loc[:, 'threat_weight'])
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'VU', 1/4, species_taxon.loc[:, 'threat_weight'])
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'NT', 1/8, species_taxon.loc[:, 'threat_weight'])
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'LC', 1/16, species_taxon.loc[:, 'threat_weight'])
  species_taxon.loc[:, 'threat_weight'] = np.where(species_taxon['category'] == 'DD', 1/4, species_taxon.loc[:, 'threat_weight'])
  
  # Get the list of species
  fnames_tr = [join(range_dir, txn, "{}.tif".format(sps.replace(" ", "_"))) 
               for sps in species_taxon['species'].tolist()]
  fnames = glob.glob(join(range_dir, txn, "*.tif"), recursive=True)
  fnames = [fname for fname in fnames if fname in fnames_tr]
  
  # Start to calculate the richness
  sr, wsr, rwsr, wrwsr = species_richness(fnames[0], species_taxon)
    
  # Add on the rest one at a time
  for fname in fnames[1:]:
    print(fname)
    _sr, _wsr, _rwsr, _wrwsr = species_richness(fname, species_taxon)
    sr = sr + _sr
    wsr = wsr + _wsr
    rwsr = rwsr + _rwsr
    wrwsr = wrwsr + _wrwsr
  
  # Get some parameters
  img = gdal.Open(fname)
  n_rows = img.RasterYSize
  n_cols = img.RasterXSize
  trans = img.GetGeoTransform()
  proj = img.GetProjection()
  ds = None
  
  # Save out
  driver = gdal.GetDriverByName("GTiff")
  
  # Species richness, sr
  sr_path = join(bio_dir, "{}_species_richness.tif".format(txn))
  out_data = driver.Create(sr_path, n_cols, n_rows, 1, gdal.GDT_Float32)
  out_data.SetGeoTransform(trans)
  out_data.SetProjection(proj)
  out_data.GetRasterBand(1).WriteArray(sr)
  out_data.FlushCache()
  out_data = None
  
  # Weighted species richness, wsr
  wsr_path = join(bio_dir, "{}_weighted_species_richness.tif".format(txn))
  out_data = driver.Create(wsr_path, n_cols, n_rows, 1, gdal.GDT_Float32)
  out_data.SetGeoTransform(trans)
  out_data.SetProjection(proj)
  out_data.GetRasterBand(1).WriteArray(wsr)
  out_data.FlushCache()
  out_data = None
  
  # Rarity-weighted species richness, rwsr
  rwsr_path = join(bio_dir, "{}_rarity_richness.tif".format(txn))
  out_data = driver.Create(rwsr_path, n_cols, n_rows, 1, gdal.GDT_Float32)
  out_data.SetGeoTransform(trans)
  out_data.SetProjection(proj)
  out_data.GetRasterBand(1).WriteArray(rwsr)
  out_data.FlushCache()
  out_data = None
  
  # Weighted rarity-weighted species richness, wrwsr
  wrwsr_path = join(bio_dir, "{}_weighted_rarity_richness.tif".format(txn))
  out_data = driver.Create(wrwsr_path, n_cols, n_rows, 1, gdal.GDT_Float32)
  out_data.SetGeoTransform(trans)
  out_data.SetProjection(proj)
  out_data.GetRasterBand(1).WriteArray(wrwsr)
  out_data.FlushCache()
  out_data = None
  
if __name__ == '__main__':
    main()
