from collections import OrderedDict
from itertools import product
import numpy as np
import re

import pandas as pd
from scipy import ndimage
import rasterio

iter_seq = snakemake.params['iterator_sequence']

params_product = product(*map(
    lambda param: snakemake.params[param],
    iter_seq
))

df = pd.DataFrame(columns=iter_seq)

landcover_categories = [
    'Crops/horticulture', 'Intensive grass', 'Extensive grass',
    'Shrubland', 'Exotic forest', 'Native forest', 'Bare ground / Others'
]
landcover_categories_ids = range(0,len(landcover_categories))

landscape_metrics_labels = ['enn_mn', 'shdi', 'contag']

# TODO parallelise?
for i, row_params in enumerate(params_product):
    landscape, topography, landclass, proportion = row_params

    data = OrderedDict({
        'landscape': landscape,
        'topography': topography,
        'landclass': landclass,
        'proportion': proportion
    })

    # Landscapes (input summary)
    re_filter = re.compile(f'{landscape}/{landscape}-landclass_{topography}_{landclass}_{proportion}')
    landscape_data = list(filter(
        re_filter.search, snakemake.params['landscapes']
    ))[0]
    with rasterio.open(landscape_data) as ds:
        band = ds.read(1, masked=True)
        total_cells = ds.height * ds.width
        num_patches = 0
        patch_sizes = []
        for class_name, id in zip(landcover_categories, landcover_categories_ids):
            per_patch_sizes = []
            # Determine percentage of total landscape occupied by class
            data[f'% {class_name}'] = np.count_nonzero(band == id) / total_cells * 100
            # Identify "patches": sub-regional areas of the same class value
            # Queen-case structure; use (2,1) for rook-case
            structure = ndimage.generate_binary_structure(2,2)
            # Label queen-case regions of landclass
            labels, num_labels = ndimage.label(np.where(band == id, 1, 0), structure=structure)
            data[f'{class_name} patches'] = num_labels
            num_patches += num_labels
            for loc in ndimage.find_objects(labels): # Slices
                # For each patch (accessed with slice), determine its size
                patch_size = np.count_nonzero(labels[loc]) * snakemake.params['resolution']**2
                patch_sizes.append(patch_size)
                per_patch_sizes.append(patch_size)
            data[f'{class_name} mean patch size'] = np.mean(per_patch_sizes)
        data['Total number of patches'] = num_patches
        data['Mean patch size'] = np.mean(patch_sizes)

    # Nitrate leaching
    re_filter = re.compile(f'/{landscape}/{topography}\.{landclass}\.{proportion}/intermediate_outputs/modified_load_n.tif')
    nitrate_modified_load_data = list(filter(
        re_filter.search, snakemake.params['ndr_intermediate_modified_load_n']
    ))[0]
    with rasterio.open(nitrate_modified_load_data) as ds:
        band = ds.read(1, masked=True)
        nLoad = band.sum()
    re_filter = re.compile(f'/{landscape}/{topography}\.{landclass}\.{proportion}/n_export.tif')
    nitrate_modified_load_data = list(filter(
        re_filter.search, snakemake.params['ndr']
    ))[0]
    with rasterio.open(nitrate_modified_load_data) as ds:
        band = ds.read(1, masked=True)
        nExport = band.sum()
    data['Total landscape N load'] = nLoad
    data['Total landscape N export'] = nExport
    data['% landscape N retained'] = (( nLoad - nExport ) / nLoad ) * 100
    data['% lanscape N exported'] = ( nExport / nLoad ) * 100

    # Soil erosion
    re_filter = re.compile(f'{landscape}/{topography}.{landclass}.{proportion}')
    erosion_data = list(filter(
        re_filter.search, snakemake.params['erosion']
    ))[0]
    with rasterio.open(erosion_data) as ds:
        band = ds.read(1, masked=True)
        data['Max per pixel soil erosion'] = band.max()
        data['Mean per pixel soil erosion'] = band.mean()
        data['Landscape sum of soil erosion'] = band.sum()/(10000/snakemake.params['resolution']**2)


    # Recreation
    recreation_data = list(filter(
        re_filter.search, snakemake.params['recreation']
    ))[0]
    with rasterio.open(recreation_data) as ds:
        band = ds.read(1, masked=True)
        data['Max per pixel recreation score'] = band.max()
        data['Mean per pixel recreation score'] = band.mean()

    # Pollination
    pollination_data = list(filter(
        re_filter.search, snakemake.params['pollinator_overlap']
    ))[0]
    with rasterio.open(pollination_data) as ds:
        band = ds.read(1, masked=True)
        total_cells = ds.height * ds.width
        data['no pollinators present'] = np.count_nonzero(band == 11) / total_cells * 100
        data['medium quality pollinator habitat present'] = np.count_nonzero(band == 16) / total_cells * 100
        data['high quality pollinator habitat present'] = np.count_nonzero(band == 21) / total_cells * 100

    # Greenhouse gases
    ghg_data = list(filter(
        re_filter.search, snakemake.params['greenhouse_gas_emissions']
    ))[0]
    with rasterio.open(ghg_data) as ds:
        band = ds.read(1, masked=True)
        data['Mean greenhouse gas emissions'] = band.mean()
        data['Landscape sum of greenhouse gas emissions'] = band.sum()/(10000/snakemake.params['resolution']**2)

    # Carbon stocks
    carbon_stock_data = list(filter(
        re_filter.search, snakemake.params['carbon_stocks']
    ))[0]
    with rasterio.open(carbon_stock_data) as ds:
        band = ds.read(1, masked=True)
        data['Mean carbon stocks'] = band.mean()
        data['Landscape sum of carbon stocks'] = band.sum()/(10000/snakemake.params['resolution']**2)

    # Landscape metrics
    landscape_metrics_data = list(filter(
        re_filter.search, snakemake.params['landscape_metrics']
    ))[0]
    ls_df = pd.read_csv(landscape_metrics_data)
    for metric in landscape_metrics_labels:
        data[metric] = ls_df[ls_df['metric'] == metric]['value'].iat[0]

    # Append row
    df = df.append(pd.Series(data=data, name=i), ignore_index=False)

df = df.reindex(columns=data.keys())
df.rename(columns={
    'landscape': 'Landscape Type',
    'topography': 'Topography ID',
    'landclass': 'Fragmentation',
    'proportion': 'Proportion set'
}, inplace=True)
df.to_excel(snakemake.output[0], sheet_name="landscape_summary")

sys.exit(0)
