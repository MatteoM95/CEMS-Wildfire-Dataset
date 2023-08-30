import os
import json
import pystac_client

import numpy as np
import pandas as pd
import planetary_computer as pc
import rasterio
import rasterio.features
import stackstac
import odc.stac

from PIL import Image

from pystac.extensions.item_assets import ItemAssetsExtension

import rasterio
from datetime import datetime

import numpy as np
import rasterio

import src.utils_variable as utils
os.environ["PLANETARY_COMPUTER_API_KEY"] = utils.planetary_computer_API_KEY

# Set the environment variable PC_SDK_SUBSCRIPTION_KEY, or set it here.
# The Hub sets PC_SDK_SUBSCRIPTION_KEY automatically.
# pc.settings.set_subscription_key(<YOUR API Key>)


class DownloadLandCover:

    def to_RGB_Mask(self, image_mask: np.ndarray, path: str, colormap: np.ndarray):

        """Save a RGB mask for the land cover

        Args:
            image_mask (np.ndarray): land cover image with only one channel
            path (str): save path
            colormap: (np.ndarray): colormap for each label 
           
        Returns:
            None
        """

        unique_values = set(np.unique(image_mask).astype(int).tolist()) # getting unique classes
        # colors = [(0,0,0), (181,254,142), (254,217,142), (254,153,41), (204,76,2), (255,255,255), (102,255,255), (204,76,2), (255,255,255), (102,255,255), (204,76,2), (255,255,255), (102,255,255)]
        colors = colormap

        Rarr = np.zeros_like(image_mask, dtype = 'uint8') # Red
        Garr = np.zeros_like(image_mask, dtype = 'uint8') # Green
        Barr = np.zeros_like(image_mask, dtype = 'uint8') # Blue
        for val, col in zip(unique_values, colors):
            Rarr[image_mask == val ] = colors[val][0] #col[0]
            Garr[image_mask == val ] = colors[val][1] #col[1]
            Barr[image_mask == val ] = colors[val][2] #col[2]

        RGB_mask = np.dstack((Rarr,Garr,Barr)) # Combining three channels
        image_RGB = Image.fromarray(RGB_mask, 'RGB')
        image_RGB.save(path)


    def download_Esri_LandCover10(self):

        """ Download and save Esri 10 class land cover mask (2020 only)
            https://planetarycomputer.microsoft.com/dataset/io-lulc#overview
            Default labels: 
            {'nodata': 0, 'water': 1, 'trees': 2, 'grass': 3, 'flooded veg': 4, 'crops': 5, 'scrub': 6, 'built area': 7, 'bare': 8, 'snow/ice': 9, 'clouds': 10}
        """

        df = pd.read_csv(utils.satelliteData)
        df['folderPath'] = utils.output_folder + df['folderPath'].astype(str)

        for index in range(len(df)):
        
            bbox = [df['left_Long'].loc[index],df['bottom_Lat'].loc[index],df['right_Long'].loc[index],df['top_Lat'].loc[index]]
            intervalDatetime = [df['interval_startDate'].loc[index], df['interval_endDate'].loc[index]]
            height_aoi = df['height'].loc[index]
            width_aoi = df['width'].loc[index]
            
            basename = df['folderPath'].loc[index].split("/")[-1]
            path_png = df['folderPath'].loc[index] + "/" + basename + "_Esri10_LC.png"
            path_tif = df['folderPath'].loc[index] + "/" + basename + "_Esri10_LC.tif"

            if os.path.exists(path_tif):
                continue

            catalog = pystac_client.Client.open(
                "https://planetarycomputer.microsoft.com/api/stac/v1",
                modifier=pc.sign_inplace,
            )

            search = catalog.search(collections=["io-lulc"], bbox=bbox, datetime=['2020-01-01T00:00:00Z', '2020-02-02T00:00:00Z'])

            items = search.item_collection()
            item = items[0]

            nodata = 0 #item.assets["data"].extra_fields["raster:bands"][0]["nodata"]

            stack = stackstac.stack(
                items, 
                epsg=4326, 
                dtype=np.ubyte, 
                fill_value=nodata, 
                bounds_latlon=bbox, 
                resolution = (abs(df['resolution_x'].loc[index]), 
                              abs(df['resolution_y'].loc[index])),
                snap_bounds=False
            )
 
            mergedRemapped = stackstac.mosaic(stack, dim="time", axis=None, nodata=0).squeeze().compute()
            merged = stackstac.mosaic(stack, dim="time", axis=None, nodata=0).squeeze().compute()

            # land cover label shifted to a lower value, negative value will be invalid data 255
            mask = mergedRemapped.data
            mask -= 1                           # shift label to lower value
            mask[mask == -1] = 255              # No data value 0 is now 255
            mergedRemapped.data = mask

            if mergedRemapped.data.shape[0] != height_aoi:
                print("Problem with height LC: " + df['folderPath'].loc[index])
            if mergedRemapped.data.shape[1] != width_aoi:
                print("Problem with width LC: " + df['folderPath'].loc[index])

            if mask.max() == 255:
                file_object = open(utils.logPath, 'a')
                file_object.write("No-data LandCover in : " + path_tif + " - " + str(mask.shape) +"\n")
                file_object.close()

            # save tif image
            mergedRemapped.rio.to_raster(raster_path=path_tif)

            class_names = merged.coords["label:classes"].item()["classes"]
            class_count = len(class_names)

            colormap = np.zeros(256)
            with rasterio.open(item.assets["data"].href) as src:
                colormap_def = src.colormap(1)  # get metadata colormap for band 1
                colormap = [ np.array(colormap_def[i]) for i in range(class_count) ]  # transform to matplotlib color format

            self.to_RGB_Mask(image_mask=merged, path = path_png, colormap = colormap )

    #
    def download_AnnualLandCover9(self):

        """ Download and save Esri 9 class annual land cover mask (2017 - 2021)
            https://planetarycomputer.microsoft.com/dataset/io-lulc-9-class
            Default labels: 
            {'No Data': 0, 'Water': 1, 'Trees': 2, 'Flooded vegetation': 4, 'Crops': 5, 'Built area': 7, 'Bare ground': 8, 'Snow/ice': 9, 'Clouds': 10, 'Rangeland': 11}
        """

        df = pd.read_csv(utils.satelliteData)
        df['folderPath'] = utils.output_folder + df['folderPath'].astype(str)

        for index in range(len(df)):
        
            bbox = [df['left_Long'].loc[index],df['bottom_Lat'].loc[index],df['right_Long'].loc[index],df['top_Lat'].loc[index]]
            intervalDatetime = [df['interval_startDate'].loc[index], df['interval_endDate'].loc[index]]
            height_aoi = df['height'].loc[index]
            width_aoi = df['width'].loc[index]

             # the land cover for 2022 is not available yet
            datetime_object = datetime.strptime(intervalDatetime[0], '%Y-%m-%dT%H:%M:%SZ')
            if datetime_object.year == 2022 or datetime_object.year == 2023:
                intervalDatetime = ['2021-11-01T00:00:00Z', '2021-12-02T00:00:00Z']

            basename = df['folderPath'].loc[index].split("/")[-1]
            path_png = df['folderPath'].loc[index] + "/" + basename + "_Annual9_LC.png"
            path_tif = df['folderPath'].loc[index] + "/" + basename + "_Annual9_LC.tif"
            
            if os.path.exists(path_tif):
                continue

            catalog = pystac_client.Client.open(
                "https://planetarycomputer.microsoft.com/api/stac/v1",
                modifier=pc.sign_inplace,
            )

            search = catalog.search(collections=["io-lulc-9-class"], bbox=bbox, datetime=intervalDatetime) #, intersects=area_of_interest)

            # Check how many items were returned
            items = search.item_collection()
            item = items[0]

            nodata = 0      #item.assets["data"].extra_fields["raster:bands"][0]["nodata"]

            stack = stackstac.stack(
                items, 
                epsg=4326, 
                dtype=np.ubyte, 
                fill_value=nodata, 
                bounds_latlon=bbox,
                resolution = (abs(df['resolution_x'].loc[index]), 
                            abs(df['resolution_y'].loc[index])),
                snap_bounds=False,
            )

            mergedRemapped = stackstac.mosaic(stack, dim="time", axis=None, nodata=0).squeeze().compute()
            merged = stackstac.mosaic(stack, dim="time", axis=None, nodata=0).squeeze().compute()

            # land cover label shifted to a lower value, negative value will be invalid data 255
            mask = mergedRemapped.data
            mask -= 1                           # shift label to lower value
            mask[mask == -1] = 255              # No data value 0 is now 255
            mergedRemapped.data = mask

            if mergedRemapped.data.shape[0] != height_aoi:
                print("Problem: " + df['folderPath'].loc[index])
            if mergedRemapped.data.shape[1] != width_aoi:
                print("Problem: " + df['folderPath'].loc[index])

            # print(df['folderPath'].loc[index] + " - " + str(mask.shape))
            if mask.max() == 255:
                file_object = open(utils.logPath, 'a')
                file_object.write("No-data LandCover in : " + path_tif + " - " + str(mask.shape) +"\n")
                file_object.close()

            # save tif image
            mergedRemapped.rio.to_raster(raster_path=path_tif)

            collection = catalog.get_collection("io-lulc-9-class")
            ia = ItemAssetsExtension.ext(collection)

            x = ia.item_assets["data"]
            class_names = {x["summary"]: x["values"][0] for x in x.properties["file:values"]}

            colormap = np.zeros(256)
            with rasterio.open(item.assets["data"].href) as src:
                colormap_def = src.colormap(1)  # get metadata colormap for band 1
                colormap = [ np.array(colormap_def[i]) for i in range(max(class_names.values())+1) ]  # transform to matplotlib color format

            self.to_RGB_Mask(image_mask=merged, path = path_png, colormap = colormap)

    
    def download_ESA_WorldCover2020(self):

        """ Download and save ESA World cover 2020 land cover mask (2020 only)
            https://planetarycomputer.microsoft.com/dataset/esa-worldcover
            Manual here: https://worldcover2020.esa.int/data/docs/WorldCover_PUM_V1.1.pdf
            Default labels: 
            {'No Data': 0, 'Tree cover ': 10, 'Shrubland': 20, 'Grassland': 30, 'Cropland': 40, 'Built-up': 50, 'Bare / sparse vegetation': 60, 'Snow/ice': 70, 'Permanent water bodies': 80, 'Herbaceous wetland': 90, 'Mangroves': 95, 'Moss and lichen': 100}
        """

        df = pd.read_csv(utils.satelliteData)
        df['folderPath'] = utils.output_folder + df['folderPath'].astype(str)
        
        for index in range(len(df)):

            bbox = [df['left_Long'].loc[index],df['bottom_Lat'].loc[index],df['right_Long'].loc[index],df['top_Lat'].loc[index]]
            intervalDatetime=[df['interval_startDate'].loc[index], df['interval_endDate'].loc[index]]
            
            basename = df['folderPath'].loc[index].split("/")[-1]
            path_png = df['folderPath'].loc[index] + "/" + basename + "_ESA_LC.png"
            path_tif = df['folderPath'].loc[index] + "/" + basename + "_ESA_LC.tif"
            
            if os.path.exists(path_tif):
                continue

            height_aoi = df['height'].loc[index]
            width_aoi = df['width'].loc[index]

            catalog = pystac_client.Client.open(
                "https://planetarycomputer.microsoft.com/api/stac/v1",
                modifier=pc.sign_inplace,
            )

            search = catalog.search(collections=["esa-worldcover"], bbox=bbox)

            # items returned from the given bbox
            items = search.item_collection()

            # UPDATE: Lines below worked before april 2023 update
            # nodata = 0      #item.assets["data"].extra_fields["raster:bands"][0]["nodata"]
            # stack = stackstac.stack(
            #     items=items, 
            #     epsg=4326, 
            #     dtype=np.ubyte, 
            #     fill_value=nodata, 
            #     bounds_latlon=bbox,
            #     resolution = (abs(df['resolution_x'].loc[index]), 
            #                   abs(df['resolution_y'].loc[index])),
            #     snap_bounds=False,
            # )
            # mergedRemapped = stackstac.mosaic(stack, dim="time", axis=None, nodata=nodata).squeeze().compute()
            # merged = stackstac.mosaic(stack, dim="time", axis=None, nodata=nodata).squeeze().compute()
           
            # UPDATE: april 2023
            # PROBLEMA: le landcover sono un pixel in altezza e uno in larghezza più grandi, usando le risoluzioni sul CSV non ci sono problemi negli altri landcover, perchè qui si???
            stack = odc.stac.load(items, 
                                crs="EPSG:4326", 
                                resolution=odc.geo.resxy_(abs(df['resolution_x'].loc[index]), \
                                                            abs(df['resolution_y'].loc[index])), 
                                bbox=bbox)

            mergedRemapped = stack["map"].isel(time=-1).load()
            merged = stack["map"].isel(time=-1).load()

            # set nodata land to 255
            mergedRemapped.attrs["nodata"] = 255
            merged.attrs["nodata"] = 255

            # reduce of some pixels (usually 1 pixel) the landcover image --> Then check if the dimensions are the same of the S2_L2A ( part of april 2023 update )
            if mergedRemapped.data.shape[0] != height_aoi:
                mergedRemapped = mergedRemapped.isel(latitude=slice(None, (height_aoi - mergedRemapped.data.shape[0])))
                merged = merged.isel(latitude=slice(None, (height_aoi - merged.data.shape[0])))
            if mergedRemapped.data.shape[1] != width_aoi:
                mergedRemapped = mergedRemapped.isel(longitude=slice(None, (width_aoi - mergedRemapped.data.shape[1])))
                merged = merged.isel(longitude=slice(None, (width_aoi - merged.data.shape[1])))

            if mergedRemapped.data.shape[0] != height_aoi:
                print(f"Problem with dimensions image height: {df['folderPath'].loc[index]} --- height_LC: {mergedRemapped.data.shape[0]} - height_L2A: {height_aoi}")
            if mergedRemapped.data.shape[1] != width_aoi:
                print(f"Problem with dimensions image width: {df['folderPath'].loc[index]} --- width_LC: {mergedRemapped.data.shape[1]} - width_L2A: {width_aoi}")

            # REMAPPING: land cover label shifted to a lower value, negative value will be invalid data 255
            mask = mergedRemapped.data
            mask[ mask == 100 ] = 110       # Moss and lichen label
            mask[ mask == 95 ] = 100        # Mangroves label
            mask = mask/10                  # single digit index
            mask -= 1                       # shift label to lower value
            mask[mask < 0 ] = 255           # No data value 0 is now 255
            mask[mask > 120] = 255          # all other values are set to nodata land 255
            mergedRemapped.data = mask

            # save tif image
            mergedRemapped.rio.to_raster(raster_path=path_tif)

            asset_href = items[0].assets["map"].href

            colormap = np.zeros(256)
            with rasterio.open(asset_href) as src:
                colormap_def = src.colormap(1)
                colormap = [ np.array(colormap_def[i]) for i in range(len(colormap_def)) ]

            self.to_RGB_Mask(image_mask=np.flip(merged.data, 0), path = path_png, colormap = colormap ) # I do not why but after april 2023 the array should be flipped