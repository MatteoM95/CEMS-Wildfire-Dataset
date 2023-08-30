import os
import numpy as np
import json

from typing import Union

import shapely
from shapely.geometry import shape, GeometryCollection
from shapely.geometry import MultiPolygon
from shapely.ops import cascaded_union, unary_union
import rasterio.features

from PIL import Image

import rasterio
import rasterio.mask
from rasterio.merge import merge

from src.cloudCover import CloudCover
import src.utils_variable as utils

# import warnings
# from shapely.errors import ShapelyDeprecationWarning
# warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)  # Avoid warning about: "__getitem__ for multi-part geometries is deprecated and will be removed in Shapely 2.0. Use the `geoms` property to access the constituent parts of a multi-part geometry.

import warnings
warnings.filterwarnings('ignore')


class RasterizeMask():

    def __init__(self)->None:
       
        self.damageLevel = utils.damageLevel
        self.copernicusPath = utils.copernicusFolderPath


    def to_RGB_tiff(self, image_tiff, path: str, contrast: int = 0, brightness: int = 0):

        if isinstance(image_tiff, rasterio.io.DatasetReader):
            image_np = image_tiff.read()
        else:
            image_np = image_tiff

        if contrast == 0:
            contrast = utils.contrast_coeff
        else:
            contrast = contrast
        
        if brightness == 0:
            brightness = utils.brightness_coeff
        else:
            brightness = brightness

        image_np = np.transpose(image_np, (1, 2, 0))
        
        # imageRGB = image_np[:,:,[3,2,1]]
        imageRGB = np.multiply(image_np[:,:,[3,2,1]], contrast) + brightness
        imageRGB[ imageRGB > 255 ] = 255
        imageRGB[ imageRGB < 0 ] = 0        
        image_RGB = Image.fromarray(imageRGB.astype(np.uint8))
        image_RGB.save(path)


    def to_RGB_Mask(self, _image_mask: np.ndarray, path: str):

        """Save a RGB mask in .png format

        Args:
            image_mask (np.ndarray): image with only one channel
            path (str): save path 
           
        Returns:
            None
        """
        image_mask = _image_mask.copy()
        image_mask[ image_mask > 200 ] = 10
        unique_values = set(np.unique(image_mask).astype(int).tolist()) # getting unique classes
        
        colors = [  (0,0,0),            # 0 = Black = No burned
                    (181,254,142),      # 1 = Greenish = Negligible damage
                    (254,217,142),      # 2 = Light Orange = Moderate damage
                    (254,153,41),       # 3 = Orange = High damage
                    (204,76,2),         # 4 = Reddish = Destruction
                    (240,20,240),       # 5 = Purple = Cloud overlap with burned area
                    (103,190,224),      # 6 = Light blue = Clear sky
                    (220,220,220),      # 7 = White = Cloud
                    (180,180,180),      # 8 = Grey = Light cloud
                    (60,60,60),         # 9 = Dark grey = Cloud's shadow
                    (255,255,255)      # 10 = White = No Data
                    ]

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


    def get_multipolygon_grading(self, EMSR_AOI: str) -> dict:

        """Read the GRA_PRODUCT_naturalLandUse.json file and create a dict with a MultiPolygon for each grading damage level

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
           
        Returns:
             dict: key: damage level, value: MultiPolygon
        """

        folderPath = self.copernicusPath + EMSR_AOI
        fileJson = ""
        multipolygon = {}
        path = ""

        for dirpath, _, filenames in os.walk(folderPath):
            for filename in [f for f in filenames if "GRA_PRODUCT_naturalLandUse" in f]:
                path = os.path.join(dirpath, filename)
                if(os.path.isfile(path)):
                    fileJson = path
                    
        if fileJson == "":
            print("File " + EMSR_AOI + " GRA_PRODUCT_naturalLandUse does not exist")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " File GRA_PRODUCT_naturalLandUse does not exist \n")
            file_object.close()
            return

        # open json GRA_PRODUCT_naturalLandUse
        with open(fileJson) as f:
                features = json.load(f)["features"]
        
        # type geometry: Multipolygon
        multipolygon_MulDestroyed = GeometryCollection([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("destruction") and ( feature["geometry"]["type"] == "Multipolygon" or feature["geometry"]["type"] == "MultiPolygon" )] ) #                                                                                                                                    p
        multipolygon_MulHigh = GeometryCollection([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("high damage") and ( feature["geometry"]["type"] == "Multipolygon" or feature["geometry"]["type"] == "MultiPolygon" ) ] )
        multipolygon_MulModerate = GeometryCollection([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("moderate damage") and ( feature["geometry"]["type"] == "Multipolygon" or feature["geometry"]["type"] == "MultiPolygon" ) ] )
        multipolygon_MulNegligible = GeometryCollection([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("negligible damage") and ( feature["geometry"]["type"] == "Multipolygon" or feature["geometry"]["type"] == "MultiPolygon" ) ] )

        # type geometry: Polygon -> Multypoligon
        multipolygons_PolDestroyed = MultiPolygon([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("destruction") and feature["geometry"]["type"] == "Polygon" ] )
        multipolygons_PolHigh = MultiPolygon([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("high damage") and feature["geometry"]["type"] == "Polygon" ] )
        multipolygons_PolModerate = MultiPolygon([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("moderate damage") and feature["geometry"]["type"] == "Polygon" ] )
        multipolygons_PolNegligible = MultiPolygon([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["properties"]["damage_gra"] in self.damageLevel.get("negligible damage") and feature["geometry"]["type"] == "Polygon" ] )

        # merge in one single Multypoligon
        multipolygon_destroyed = unary_union([multipolygons_PolDestroyed, multipolygon_MulDestroyed])
        if isinstance(multipolygon_destroyed, shapely.geometry.polygon.Polygon):
            multipolygon_destroyed = multipolygons_PolDestroyed

        multipolygon_high = unary_union([multipolygons_PolHigh, multipolygon_MulHigh])
        if isinstance(multipolygon_high, shapely.geometry.polygon.Polygon):
            multipolygon_high = multipolygons_PolHigh

        multipolygon_moderate = unary_union([multipolygons_PolModerate, multipolygon_MulModerate])
        if isinstance(multipolygon_moderate, shapely.geometry.polygon.Polygon):
            multipolygon_moderate = multipolygons_PolModerate

        multipolygon_negligible = unary_union([multipolygons_PolNegligible, multipolygon_MulNegligible])
        if isinstance(multipolygon_negligible, shapely.geometry.polygon.Polygon):
            multipolygon_negligible = multipolygons_PolNegligible

        # Check if the object is now a Multypoligon
        if isinstance(multipolygon_negligible, shapely.geometry.polygon.Polygon):
            print(EMSR_AOI + " TypeError: object of type 'Polygon' has no len() ")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " multipolygon_negligible -> TypeError: object of type 'Polygon' has no len() \n")
            file_object.close()
            multipolygon_negligible = []
        if isinstance(multipolygon_moderate, shapely.geometry.polygon.Polygon):
            print(EMSR_AOI + " TypeError: object of type 'Polygon' has no len() ")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " multipolygon_moderate -> TypeError: object of type 'Polygon' has no len() \n")
            file_object.close()
            multipolygon_moderate = []
        if isinstance(multipolygon_high, shapely.geometry.polygon.Polygon):
            print(EMSR_AOI + " TypeError: object of type 'Polygon' has no len() ")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " multipolygon_high -> TypeError: object of type 'Polygon' has no len() \n")
            file_object.close()
            multipolygon_high = []
        if isinstance(multipolygon_destroyed, shapely.geometry.polygon.Polygon):
            print(EMSR_AOI + " TypeError: object of type 'Polygon' has no len() ")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " multipolygon_destroyed -> TypeError: object of type 'Polygon' has no len() \n")
            file_object.close()
            multipolygon_destroyed = []

        # create the dict with different level of damage
        if len(multipolygon_negligible) != 0:
            multipolygon["negligible"] = multipolygon_negligible
        if len(multipolygon_moderate) != 0:
            multipolygon["moderate"] = multipolygon_moderate
        if len(multipolygon_high) != 0:
            multipolygon["high"] = multipolygon_high
        if len(multipolygon_destroyed) != 0:
            multipolygon["destroyed"] = multipolygon_destroyed
        
        if multipolygon is None:
            print(EMSR_AOI + " GRA TypeError: object of type 'NoneType'  ")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " GRA multipolygon -> TypeError: object of type 'NoneType' \n")
            file_object.close()
            return None

        if len(multipolygon) == 0:
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " - GRA Multipolygon len = 0 \n")
            file_object.close()
            return None

        return multipolygon

    
    def get_multipolygon_observedEvent(self, EMSR_AOI: str, observedEvent: str):

        """Read the FEP_PRODUCT_observedEventA.json or DEL_PRODUCT_observedEventA.json file and create a MultiPolygon for the burned area

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            observedEvent (str): FEP or DEL PRODUCT_observedEventA 
           
        Returns:
             MultiPolygon: geometry that defines the burned area
        """

        folderPath = self.copernicusPath + EMSR_AOI
        fileJson = ""
        path = ""

        for dirpath, _, filenames in os.walk(folderPath):
            for filename in [f for f in filenames if observedEvent in f]:
                path = os.path.join(dirpath, filename)
                if(os.path.isfile(path)):
                    fileJson = path

        if fileJson == "":
            print("File " + EMSR_AOI + " " + observedEvent + " does not exist")
            file_object = open(utils.logPath, 'a')
            file_object.write(EMSR_AOI + " File " + observedEvent + "  does not exist\n")
            file_object.close()
            return

        with open(fileJson) as f:
                features = json.load(f)["features"]

        # type geometry: Multipolygon and Polygon
        multipolygon_MulDEL = GeometryCollection([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and ( feature["geometry"]["type"] == "Multipolygon" or feature["geometry"]["type"] == "MultiPolygon" ) ] )
        multipolygon_PolDEL = MultiPolygon([shape(feature["geometry"]).buffer(0)  for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon" ] )

        # merge in one single Multypoligon
        if len(multipolygon_MulDEL) == 0 and len(multipolygon_PolDEL) == 0:
            return
        elif len(multipolygon_MulDEL) == 0 and len(multipolygon_PolDEL) != 0:
            multipolygon_burned = multipolygon_PolDEL
        elif len(multipolygon_MulDEL) != 0 and len(multipolygon_PolDEL) == 0:
            multipolygon_burned = multipolygon_MulDEL
        else:
            multipolygon_burned = cascaded_union([multipolygon_MulDEL, multipolygon_PolDEL])
        
        return multipolygon_burned
    

    def get_Tiff_Merged(self, output_folder: str, EMSR_AOI: str, files: list):
        
        """Merge in a single tiff image a list of images

        Args:
            files (list): list of paths to tiff images
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder 
           
        Returns:
            np.ndarray: merged tiff image
        """

        tif_list = []

        for num in range(len(files)):                
            src = rasterio.open(output_folder + EMSR_AOI + files[num])
            tif_list.append(src)

        # merge in one single tiff
        tif_merged, out_trans = merge(tif_list)
        out_meta = src.meta.copy()

        out_meta.update(
                    {
                        "height": tif_merged.shape[1],
                        "width": tif_merged.shape[2],
                        "transform": out_trans,
                    }
                )

        if utils.saveMergedTif:
            with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + '_merged.tif', "w", **out_meta) as dest:
                    dest.write(tif_merged)

        self.to_RGB_tiff(tif_merged, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + '_merged.png')

        return tif_merged
    
    
    def rasterizeMaskGrading(self, output_folder: str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry create the related raster image for grading mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """
        
        multipolygon = []
        multipolygon = self.get_multipolygon_grading(EMSR_AOI)

        if multipolygon is None:
            return

        for num in range(len(files)):
            
            path = output_folder + EMSR_AOI + files[num].split('/')[0] + "/"
            
            with rasterio.open(output_folder + EMSR_AOI + files[num]) as tif_image:
                
                self.to_RGB_tiff(tif_image, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + '_S2L2A.png')

                # cloud cover mask
                _, cloudMask = self.createCloudMask(tif_image, 
                                                    path = path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2),
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = False) 

                # creation grading mask 
                if "destroyed" in multipolygon:
                    out_image_destroyed, out_transform = rasterio.mask.mask(tif_image, multipolygon.get("destroyed"))
                    mask_destroyed = out_image_destroyed[2,:,:]+out_image_destroyed[3,:,:]+out_image_destroyed[1,:,:]>0
                    out_image_destroyed_masked = np.where(~mask_destroyed, np.zeros(out_image_destroyed[0].shape), 4) #severity 4
                    if utils.saveGradingRGB and utils.saveRGBImage:              
                        self.to_RGB_Mask(out_image_destroyed_masked, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1]  + "_" + str(num+1).zfill(2) + "_GRA_DamageMaskDestroyed.png")
                else:
                    out_image_destroyed_masked = 0

                if "high" in multipolygon:
                    out_image_high, out_transform = rasterio.mask.mask(tif_image, multipolygon.get("high"))
                    mask_high = out_image_high[2,:,:]+out_image_high[3,:,:]+out_image_high[1,:,:]>0
                    out_image_high_masked = np.where(~mask_high, np.zeros(out_image_high[0].shape), 3) #severity 3
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_high_masked, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1]  + "_" + str(num+1).zfill(2) + "_GRA_DamageMaskHigh.png")
                else:
                    out_image_high_masked = 0

                if "moderate" in multipolygon:
                    out_image_moderate, out_transform = rasterio.mask.mask(tif_image, multipolygon.get("moderate"))
                    mask_moderate = out_image_moderate[2,:,:]+out_image_moderate[3,:,:]+out_image_moderate[1,:,:]>0
                    out_image_moderate_masked = np.where(~mask_moderate, np.zeros(out_image_moderate[0].shape), 2) #severity 2
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_moderate_masked, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1]  + "_" + str(num+1).zfill(2) + "_GRA_DamageMaskModerate.png")
                else:
                    out_image_moderate_masked = 0

                if "negligible" in multipolygon:
                    out_image_negligible, out_transform = rasterio.mask.mask(tif_image, multipolygon.get("negligible"))
                    mask_negligible = out_image_negligible[2,:,:]+out_image_negligible[3,:,:]+out_image_negligible[1,:,:]>0
                    out_image_negligible_masked = np.where(~mask_negligible, np.zeros(out_image_negligible[0].shape), 1) #severity 1
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_negligible_masked, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA_DamageMaskNegligible.png")
                else:
                    out_image_negligible_masked = 0
                
                # merge in one single mask all levels
                final_mask = out_image_high_masked + out_image_moderate_masked + out_image_negligible_masked + out_image_destroyed_masked

                # in case of overlap between masks, the damage level is set to moderate damage 
                final_mask[final_mask >= 5] = 3

                # save mask in ASCII
                # np.savetxt(path + "/mask.txt", final_mask.astype(int), fmt='%i') 

                # validity mask
                image = tif_image.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                final_mask[validity_mask == 255] = 255

                # save grading as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(final_mask, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA.png")

                out_meta = tif_image.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8
                    } 
                )

                with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(final_mask, axis=0) )

                # Make delineation mask if not exist
                if not os.path.exists(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif"):
                    delineation_mask = np.where(final_mask > 0, 1, 0)
                    delineation_mask[validity_mask == 255] = 255
                    if utils.saveRGBImage:
                        self.to_RGB_Mask(delineation_mask*3, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.png")
                    with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(delineation_mask, axis=0))

                # if delineation mask exists, save the largest one
                if os.path.exists(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif"):
                    with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif") as src:
                        delineation_mask = src.read()
                        delineation_mask = np.squeeze(delineation_mask)
                        countBurnedDEL = np.count_nonzero(delineation_mask > 0)
                        countBurnedGRA = np.count_nonzero(final_mask > 0)
                        if countBurnedDEL < countBurnedGRA:
                            delineation_mask = np.where(final_mask > 0, 1, 0)
                            delineation_mask[validity_mask == 255] = 255
                            with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif", "w", **out_meta) as dest:
                                dest.write(np.expand_dims(delineation_mask, axis=0))
                            if utils.saveRGBImage:
                                self.to_RGB_Mask(delineation_mask*3, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.png")


                # save grading image cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( final_mask, cloudMask )
                    
                    with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA_CleanDamageMask.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(final_mask, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA_OverlapMask.png")
                        self.to_RGB_Mask( np.multiply(final_mask, cleanMask) , path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_GRA_CleanDamageMask.png")
    

    def rasterizeMaskDelineation(self, output_folder: str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry create the related raster image for delineation mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """

        multipolygon_burned = self.get_multipolygon_observedEvent(EMSR_AOI, "DEL_PRODUCT_observedEventA")

        for num in range(len(files)):

            path = output_folder + EMSR_AOI + files[num].split('/')[0] + "/"

            with rasterio.open(output_folder + EMSR_AOI + files[num]) as tif_image:

                self.to_RGB_tiff(tif_image, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + '_S2L2A.png')

                #cloud cover mask
                _, cloudMask = self.createCloudMask(tif_image, 
                                                    path = path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2),
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = False)

                out_image_burned, out_transform = rasterio.mask.mask(tif_image, multipolygon_burned)
                mask_burned = out_image_burned[2,:,:]+out_image_burned[3,:,:]+out_image_burned[1,:,:]>0
                out_image_burned_masked = np.where(~mask_burned, np.zeros(out_image_burned[0].shape), 1) # level burned area set to 1
                
                # validity mask
                image = tif_image.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                out_image_burned_masked[validity_mask == 255] = 255

                # save delineation as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(out_image_burned_masked*3, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.png")

                out_meta = tif_image.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )

                with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(out_image_burned_masked, axis=0))
                
                # save delineation image cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( out_image_burned_masked, cloudMask )

                    with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL_CleanDelineationMask.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(out_image_burned_masked, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL_OverlapMask.png")
                        self.to_RGB_Mask( np.multiply(out_image_burned_masked, cleanMask) , path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_DEL_CleanDelineationMask.png")


    def rasterizeMaskFEP(self, output_folder: str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry create the related raster image for first estimation mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """

        multipolygon_burned = self.get_multipolygon_observedEvent(EMSR_AOI, "FEP_PRODUCT_observedEventA")

        for num in range(len(files)):

            path = output_folder + EMSR_AOI + files[num].split('/')[0] + "/"

            with rasterio.open(output_folder + EMSR_AOI + files[num]) as tif_image:
                
                self.to_RGB_tiff(tif_image, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + '_S2L2A.png')

                #cloud cover mask
                _, cloudMask = self.createCloudMask(tif_image, 
                                                    path = path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2),
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = False)

                out_image_burned, out_transform = rasterio.mask.mask(tif_image, multipolygon_burned)
                mask_burned= out_image_burned[2,:,:]+out_image_burned[3,:,:]+out_image_burned[1,:,:]>0
                out_image_burned_masked = np.where(~mask_burned, np.zeros(out_image_burned[0].shape), 1) # level burned area set to 1
                
                # validity mask
                image = tif_image.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                out_image_burned_masked[validity_mask == 255] = 255

                # save first estimation as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(out_image_burned_masked*3, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_FEP.png")

                out_meta = tif_image.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )

                with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_FEP.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(out_image_burned_masked, axis=0))
                
                # save first estimation image cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( out_image_burned_masked, cloudMask )

                    with rasterio.open(path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_FEP_CleanEstimationMask.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(out_image_burned_masked, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_FEP_OverlapMask.png")
                        self.to_RGB_Mask( np.multiply(out_image_burned_masked, cleanMask) , path + EMSR_AOI.split('/')[0] + "_" + EMSR_AOI.split('/')[1] + "_" + str(num+1).zfill(2) + "_FEP_CleanEstimationMask.png")


    def rasterizeMaskGradingMerged(self, output_folder: str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry merge in one single image tiff the list of images and save the related raster image for grading mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """

        multipolygon = []
        multipolygon = self.get_multipolygon_grading(EMSR_AOI)

        if multipolygon is None:
            return

        tif_merged = self.get_Tiff_Merged(output_folder, EMSR_AOI, files)

        if utils.saveMergedTif:
            with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + '_merged.tif') as tif_merged:

                #cloud cover mask
                _, cloudMask = self.createCloudMask(tif_merged, 
                                                    path = output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1],
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = True)

                if "destroyed" in multipolygon:
                    out_image_destroyed, out_transform = rasterio.mask.mask(tif_merged, multipolygon.get("destroyed"))
                    mask_destroyed = out_image_destroyed[2,:,:]+out_image_destroyed[3,:,:]+out_image_destroyed[1,:,:]>0
                    out_image_destroyed_masked = np.where(~mask_destroyed, np.zeros(out_image_destroyed[0].shape), 4) # severity 4
                    if utils.saveGradingRGB and utils.saveRGBImage:               
                        self.to_RGB_Mask(out_image_destroyed_masked, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_DamageMaskDestroyed_merged.png")
                else:
                    out_image_destroyed_masked = 0

                if "high" in multipolygon:
                    out_image_high, out_transform = rasterio.mask.mask(tif_merged, multipolygon.get("high"))
                    mask_high = out_image_high[2,:,:]+out_image_high[3,:,:]+out_image_high[1,:,:]>0
                    out_image_high_masked = np.where(~mask_high, np.zeros(out_image_high[0].shape), 3) # severity 3
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_high_masked, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_DamageMaskHigh_merged.png")
                else:
                    out_image_high_masked = 0

                if "moderate" in multipolygon:
                    out_image_moderate, out_transform = rasterio.mask.mask(tif_merged, multipolygon.get("moderate"))
                    mask_moderate = out_image_moderate[2,:,:]+out_image_moderate[3,:,:]+out_image_moderate[1,:,:]>0
                    out_image_moderate_masked = np.where(~mask_moderate, np.zeros(out_image_moderate[0].shape), 2) # severity 2
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_moderate_masked, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_DamageMaskModerate_merged.png")
                else:
                    out_image_moderate_masked = 0

                if "negligible" in multipolygon:
                    out_image_negligible, out_transform = rasterio.mask.mask(tif_merged, multipolygon.get("negligible"))
                    mask_negligible = out_image_negligible[2,:,:]+out_image_negligible[3,:,:]+out_image_negligible[1,:,:]>0
                    out_image_negligible_masked = np.where(~mask_negligible, np.zeros(out_image_negligible[0].shape), 1) # severity 1
                    if utils.saveGradingRGB and utils.saveRGBImage:
                        self.to_RGB_Mask(out_image_negligible_masked, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_DamageMaskNegligible_merged.png")
                else:
                    out_image_negligible_masked = 0

                out_meta = tif_merged.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )
                
                # merge in one single mask all levels
                final_mask = out_image_high_masked + out_image_moderate_masked + out_image_negligible_masked + out_image_destroyed_masked
                
                # in case of overlap between mask, the damage level is set to moderate damage 
                final_mask[final_mask >= 5] = 3

                # save mask in ASCII
                # np.savetxt(ESMR_dataPath + "mask.txt", final_mask.astype(int), fmt='%i') 
                
                # validity mask
                image = tif_merged.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                final_mask[validity_mask == 255] = 255

                # save grading as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(final_mask, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_merged.png")
                
                with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_merged.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(final_mask, axis=0))
                
                # save grading mask cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( final_mask, cloudMask )
                
                    with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_CleanDamageMask_merged.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(final_mask, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_OverlapMask_merged.png")
                        self.to_RGB_Mask( np.multiply(final_mask, cleanMask) , output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_GRA_CleanDamageMask_merged.png")


    def rasterizeMaskDelineationMerged(self, output_folder:str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry merge in one single image tiff the list of images and save the related raster image for delineation mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """

        multipolygon_burned = self.get_multipolygon_observedEvent(EMSR_AOI, "DEL_PRODUCT_observedEventA")
             
        tif_merged = self.get_Tiff_Merged(output_folder, EMSR_AOI, files)

        if utils.saveMergedTif:
            with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + '_merged.tif') as tif_merged:

                # cloud cover mask
                _, cloudMask = self.createCloudMask(tif_merged, 
                                                    path = output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1],
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = True) 

                out_image_burned, out_transform = rasterio.mask.mask(tif_merged, multipolygon_burned)
                mask_burned = out_image_burned[2,:,:] + out_image_burned[3,:,:] + out_image_burned[1,:,:] > 0
                out_image_burned_masked = np.where(~mask_burned, np.zeros(out_image_burned[0].shape), 1) # level burned area set to 1

                # validity mask
                image = tif_merged.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                out_image_burned_masked[validity_mask == 255] = 255

                # save delineation as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(out_image_burned_masked*3, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_DEL_merged.png")

                out_meta = tif_merged.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )

                with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_DEL_merged.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(out_image_burned_masked, axis=0))
                
                # save delineation mask cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( out_image_burned_masked, cloudMask )
                
                    with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_DEL_CleanDelineationMask_merged.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(out_image_burned_masked, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_DEL_overlapMask_merged.png")
                        self.to_RGB_Mask( np.multiply(out_image_burned_masked, cleanMask) , output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_DEL_CleanDelineationMask_merged.png")


    def rasterizeMaskFEPMerged(self, output_folder:str, EMSR_AOI: str, files: list):

        """From the MultiPolygon geometry merge in one single image tiff the list of images and save the related raster image for first estimation mask

        Args:
            EMSR_AOI (str): EMSR activation and Area of Interest
            output_folder (str): output folder
            files (list): list of paths to tiff images
        """

        multipolygon_burned = self.get_multipolygon_observedEvent(EMSR_AOI, "FEP_PRODUCT_observedEventA")
            
        tif_merged = self.get_Tiff_Merged(output_folder, EMSR_AOI, files)

        if utils.saveMergedTif:
            with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + '_merged.tif') as tif_merged:

                # cloud cover mask
                _, cloudMask = self.createCloudMask(tif_merged, 
                                                    path = output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1],
                                                    saveCloudMaskToFile = True, 
                                                    saveCloudTiffToFile = True, 
                                                    merged = True)

                out_image_burned, out_transform = rasterio.mask.mask(tif_merged, multipolygon_burned)
                mask_burned= out_image_burned[2,:,:]+out_image_burned[3,:,:]+out_image_burned[1,:,:]>0
                out_image_burned_masked = np.where(~mask_burned, np.zeros(out_image_burned[0].shape), 1) # level burned area set to 1

                # validity mask
                image = tif_merged.read()
                validity_mask = np.where(np.all(image == 0, axis=0), 255, 0)

                out_image_burned_masked[validity_mask == 255] = 255

                # save first estimate as png image
                if utils.saveRGBImage:
                    self.to_RGB_Mask(out_image_burned_masked*3, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_FEP_merged.png")

                out_meta = tif_merged.meta

                out_meta.update(
                    {
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )

                with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_FEP_merged.tif", "w", **out_meta) as dest:
                    dest.write(np.expand_dims(out_image_burned_masked, axis=0))
                
                # save first estimation mask cleans from clouds and the overlap between clouds and bruned area
                if utils.saveOverlapMask:
                    _, _, _, _, _, overlapMask, cleanMask = self.statisticsOverPixel( out_image_burned_masked, cloudMask )
                
                    with rasterio.open(output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_FEP_CleanEstimationMask_merged.tif", "w", **out_meta) as dest:
                        dest.write(np.expand_dims(np.multiply(out_image_burned_masked, cleanMask), axis=0))

                    if utils.saveRGBImage:
                        self.to_RGB_Mask( overlapMask*5, output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_FEP_overlapMask_merged.png")
                        self.to_RGB_Mask( np.multiply(out_image_burned_masked, cleanMask) , output_folder + EMSR_AOI + EMSR_AOI.split("/")[0] + "_" + EMSR_AOI.split("/")[1] + "_FEP_CleanEstimationMask_merged.png")

    
    def createCloudMask(self, 
                        image: Union[rasterio.io.DatasetReader, str], 
                        path: str = "",
                        saveCloudMaskToFile: bool = False, 
                        saveCloudTiffToFile: bool = False, 
                        merged: bool = False):
        
        """Create the cloud mask using cloudSen12's model and save it 

        Args:
            image (rasterio.io.DatasetReader): image tiff where predict cloud coverage
            path (str): where to save the cloud mask
            saveCloudMaskToFile (bool): if True save the cloud mask 
            saveCloudTiffToFile (bool): if True save the cloud prediction with 4 layers
            merged (bool): if the image passed as argument is merged
           
        Returns:
             np.ndarray, np.ndarray: return the prediction from cloudSen12 and the cloud mask
        """

        if isinstance(image, str):
            return self.createCloudMaskfromFilePath(image, path, saveCloudMaskToFile, saveCloudTiffToFile, merged) 

        cloudCov = CloudCover(cropsize = 512, overlap = 0)
        resultMask = cloudCov.predictCover(image.read()*10000)

        if utils.customCloudPercentage:
            cloudMaskCloud = resultMask[1] > utils.cloudPercentage 
            cloudMaskLightCloud = resultMask[2] > utils.lightCloudPercentage
            cloudMask = cloudMaskCloud | cloudMaskLightCloud #( cloudMaskCloud > 0 ) & ( cloudMaskLightCloud > 0 )
        else:
            cloudMask = np.argmax(resultMask.to_numpy(), axis=0)
            
        if merged:
            pathMask = path + "_CM_merged.tif"
            pathCloud = path + "_Cloud_merged.tif"
            pathRGB = path + "_CM_merged.png"
        else:
            pathMask = path + "_CM.tif"
            pathCloud = path + "_Cloud.tif"
            pathRGB = path + "_CM.png"

        if saveCloudMaskToFile:
            out_metaMask = image.meta

            out_metaMask.update({   
                    "count": 1, # save 1 level for mask
                    "dtype": rasterio.uint8,
                }
            )
 
            with rasterio.open(pathMask, "w", **out_metaMask) as dest:
                dest.write(np.expand_dims(cloudMask, axis=0))
            
        if saveCloudTiffToFile:
            out_metaCloud = image.meta

            out_metaCloud.update({   
                    "count": 4, # save 4 level for clouds
                    "dtype": rasterio.uint16
                }
            )
            
            if utils.saveCloudImage:
                with rasterio.open(pathCloud, "w", **out_metaCloud) as dest:
                    dest.write(resultMask)
            
        if utils.saveRGBImage:
            self.to_RGB_Mask(cloudMask + 6, pathRGB)

        return resultMask, cloudMask
    

    def createCloudMaskfromFilePath(self,
                                    image: str, 
                                    path: str,
                                    saveCloudMaskToFile: bool = True,
                                    saveCloudTiffToFile: bool = False,
                                    ):
        """Create the cloud mask using cloudSen12's model from a given file path and save it 

        Args:
            image (str): name of image tiff where predict cloud coverage
            path (str): where to save the cloud mask
            saveCloudMaskToFile (bool): if True save the cloud mask 
            saveCloudTiffToFile (bool): if True save the cloud prediction with 4 layers
           
        Returns:
             np.ndarray, np.ndarray: return the prediction from cloudSen12 and the cloud mask
        """

        if not os.path.exists(path):
            print(f"File path wrong. File {path} does not exist")
            return 0, 0

        if not path.endswith("_S2L2A.tiff"):
            path = path + "/" + image

        with rasterio.open(path) as tiff:

            tif_image = tiff.read()

            cloudCov = CloudCover(cropsize = 512, overlap = 0)

            resultMask = cloudCov.predictCover(tif_image*10000)

            # layer 0: clean sky / sunny
            # layer 1: cloudy / heavy clouds
            # layer 2: thin layer / fog / smoke / light clouds
            # layer 3: cloud's shadow

            if utils.customCloudPercentage:
                cloudMaskCloud = resultMask[1] > utils.cloudPercentage 
                cloudMaskLightCloud = resultMask[2] > utils.lightCloudPercentage
                cloudMask = cloudMaskCloud | cloudMaskLightCloud
            else:
                cloudMask = np.argmax(resultMask.to_numpy(), axis=0)

            pathMask = path.replace("_S2L2A.tiff", "_CM.tif")
            pathCloud = path.replace("_S2L2A.tiff", "_Cloud.tif")
            pathRGB = path.replace("_S2L2A.tiff", "_CM.png")

            if saveCloudMaskToFile:
                out_metaMask = tiff.meta

                out_metaMask.update({   
                        "count": 1, # save 1 level for mask
                        "dtype": rasterio.uint8,
                    }
                )

                with rasterio.open(pathMask, "w", **out_metaMask) as dest:
                    dest.write(np.expand_dims(cloudMask, axis=0))
            
            # save prediction from cloudSen12
            if utils.saveCloudImage or saveCloudTiffToFile:
                out_metaCloud = tiff.meta

                out_metaCloud.update({   
                        "count": 4, # save 4 levels for clouds
                        "dtype": rasterio.uint16
                    }
                )
                
                with rasterio.open(pathCloud, "w", **out_metaCloud) as dest:
                    dest.write(resultMask)

            # save cloud mask with different levels as png image
            if utils.saveRGBImage:
                self.to_RGB_Mask(cloudMask + 6, pathRGB)
            
            return resultMask, cloudMask


    def writeCloudcoverageCSV(self):

        """Iterate over all dataset and retrieve the information about cloud coverage and save in a .csv"""

        print("Writing cloudCoverage.csv....")

        for dirpath, _, filenames in os.walk(utils.output_folder):
            for filename in [f for f in filenames if "_CM.tif" in f]:
                
                pathToCloudMask = dirpath + "/" + filename
                pathToBurnedGRA = dirpath + "/" + filename.replace("_CM", "_GRA")
                pathToBurnedDEL = dirpath + "/" + filename.replace("_CM", "_DEL")
                pathToBurnedFEP = dirpath + "/" + filename.replace("_CM", "_FEP")
                pathToJson = dirpath + "/" + filename.replace("_CM.tif", "_S2L2A.json")
                
                Emsr_AOI = dirpath.split("/")[-1]

                with open(pathToJson) as f:
                    json_file = json.load(f)
                    payload = json_file["payload"]
                    startDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["from"]
                    endDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["to"]

                if os.stat(utils.cloudCoveragecsv).st_size == 0:
                    file_object = open(utils.cloudCoveragecsv, 'a')
                    file_object.write("EMSR_AOI,folderPath,startDate,endDate,height,width,sizeImage,burnedPixel,cloudPixel,countOverlap,percentageCloud,percentageOverlap,Type\n")
                    file_object.close()

                if os.path.exists(pathToBurnedGRA):
                    cloudMaskDatasetReader = rasterio.open(pathToCloudMask)
                    burnedMasDatasetReader = rasterio.open(pathToBurnedGRA)
                    cloudMask = np.squeeze(cloudMaskDatasetReader.read())
                    burnedMask = np.squeeze(burnedMasDatasetReader.read())
                    countBurned, countCloud, countOverlapped, percentageCloud, percentageOverlap, _, _ = self.statisticsOverPixel( burnedMask = burnedMask, cloudMask = cloudMask)
                    file_object = open(utils.cloudCoveragecsv, 'a')
                    file_object.write(Emsr_AOI + "," + dirpath + "," + str(startDate) + "," + str(endDate) + "," + str(burnedMask.shape[0]) + "," + str(burnedMask.shape[1]) + "," + str(burnedMask.shape[1] * burnedMask.shape[0]) + "," + str(round(countBurned, 2)) + "," + str(round(countCloud, 2)) + "," + str(round(countOverlapped, 2)) + "," + str(round(percentageCloud, 2)) + "," + str(round(percentageOverlap, 2)) + ",GRA\n")
                    file_object.close()
                    cloudMaskDatasetReader.close()
                    burnedMasDatasetReader.close()

                if os.path.exists(pathToBurnedDEL):
                    cloudMaskDatasetReader = rasterio.open(pathToCloudMask)
                    burnedMasDatasetReader = rasterio.open(pathToBurnedDEL)
                    cloudMask = np.squeeze(cloudMaskDatasetReader.read())
                    burnedMask = np.squeeze(burnedMasDatasetReader.read())
                    countBurned, countCloud, countOverlapped, percentageCloud, percentageOverlap, _, _ = self.statisticsOverPixel( burnedMask = burnedMask, cloudMask = cloudMask)
                    file_object = open(utils.cloudCoveragecsv, 'a')
                    file_object.write(Emsr_AOI + "," + dirpath + "," + str(startDate) + "," + str(endDate) + "," + str(burnedMask.shape[0]) + "," + str(burnedMask.shape[1]) + "," + str(burnedMask.shape[1] * burnedMask.shape[0]) + "," + str(round(countBurned, 2)) + "," + str(round(countCloud, 2)) + "," + str(round(countOverlapped, 2)) + "," + str(round(percentageCloud, 2)) + "," + str(round(percentageOverlap, 2)) + ",DEL\n")
                    file_object.close()
                    cloudMaskDatasetReader.close()
                    burnedMasDatasetReader.close()

                if os.path.exists(pathToBurnedFEP):
                    cloudMaskDatasetReader = rasterio.open(pathToCloudMask)
                    burnedMasDatasetReader = rasterio.open(pathToBurnedFEP)
                    cloudMask = np.squeeze(cloudMaskDatasetReader.read())
                    burnedMask = np.squeeze(burnedMasDatasetReader.read())
                    countBurned, countCloud, countOverlapped, percentageCloud, percentageOverlap, _, _ = self.statisticsOverPixel( burnedMask = burnedMask, cloudMask = cloudMask)
                    file_object = open(utils.cloudCoveragecsv, 'a')
                    file_object.write(Emsr_AOI + "," + dirpath + "," + str(startDate) + "," + str(endDate) + "," + str(burnedMask.shape[0]) + "," + str(burnedMask.shape[1]) + "," + str(burnedMask.shape[1] * burnedMask.shape[0]) + "," + str(round(countBurned, 2)) + "," + str(round(countCloud, 2)) + "," + str(round(countOverlapped, 2)) + "," + str(round(percentageCloud, 2)) + "," + str(round(percentageOverlap, 2)) + ",FEP\n")
                    file_object.close()
                    cloudMaskDatasetReader.close()
                    burnedMasDatasetReader.close()


    def statisticsOverPixel(self, burnedMask: np.ndarray, cloudMask: np.ndarray):

        """Calculate the statistics about burned area and cloud coverage.   

        Args:
            burnedMask (np.ndarray): burned mask ( shape: [H, W] )
            cloudMask (np.ndarray): cloud mask ( shape: [H, W] )
           
        Returns:
            countBurned (int): number of burned pixel in the mask
            countCloud (int): number of cloud pixel in the mask
            countOverlapped (int): number of pixel in the intersection between cloud and burned
            percentageCloud (float): percentage of the mask covered by the cloud
            percentageOverlap (float): percentage of the mask with overlapping pixels between burned and cloud mask 
            overlapMask (np.ndarray): intresection mask between burned mask and cloud mask
            cleanMask (np.ndarray): burned mask with removed values covered by cloud mask
        """

        countCloud = 0
        countOverlapped = 0
        
        if utils.customCloudPercentage:
            overlapMask = ( cloudMask > 0 ) & ( burnedMask > 0 ) # intersection between cloud mask and burned area mask
            cleanMask = ( burnedMask > 0 ) ^ ( ( cloudMask > 0 ) & ( burnedMask > 0 ) ) # area burned clean from clouds

            countBurned = np.count_nonzero(burnedMask > 0)
            countCloud = np.count_nonzero(cloudMask > 0)
            countOverlapped = np.count_nonzero(overlapMask > 0)
        else:
            overlapMask = ( cloudMask == 1 ) & ( burnedMask > 0 ) # intersection between cloud mask and burned area mask
            cleanMask = ( burnedMask > 0 ) ^ ( ( cloudMask == 1 ) & ( burnedMask > 0 ) ) # area burned clean from clouds

            countBurned = np.count_nonzero(burnedMask > 0)
            countCloud = np.count_nonzero(cloudMask == 1)
            countOverlapped = np.count_nonzero(overlapMask > 0)
            # print(f"ShapeBurned: {burnedMask.shape}, shapeCloud: {cloudMask.shape}, burned: {countBurned}, cloud: {countCloud}, overlap: {countOverlapped}")

        # if len(burnedMask.shape) > 2:
        #     burnedMask = np.transpose(burnedMask, (1, 2, 0)) # channel last

        percentageCloud = countCloud/( burnedMask.shape[0]*burnedMask.shape[1] ) * 100

        if countBurned != 0:
            percentageOverlap = countOverlapped/countBurned * 100
        else: 
            percentageOverlap = 0

        return countBurned, countCloud, countOverlapped, percentageCloud, percentageOverlap, overlapMask, cleanMask
