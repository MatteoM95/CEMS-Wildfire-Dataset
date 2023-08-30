import pandas as pd
import numpy as np
import json
import os
import shutil

from datetime import datetime, timedelta
from decimal import *
import rasterio.features
from geopy.geocoders import Nominatim

import shapely
from shapely.geometry import shape, GeometryCollection, MultiPolygon

import warnings
warnings.filterwarnings('ignore')

import src.utils_variable as utils
from src.retrivial import ImageRetriever
from src.rasterizeMask import RasterizeMask
from src.downloadLandCover import DownloadLandCover


class DownloadCEMSDataset:
    
    def __init__(self):

        with open(utils.copernicusFolderPath + utils.dateJson) as json_file:
            self.dataDict = json.load(json_file)


    def download_EMSR_Manager(self, 
                              allActivation: bool = False, 
                              Emsr: str = None,
                              Aoi: str = None,
                              grading:bool = True, 
                              delineation:bool = True, 
                              estimation: bool = True):
        
        """Manage the download and creation of the dataset CEMS and which activation convert to raster, 
        create cloud cover and donwload the relative land cover.

        Args:
            allActivation (bool): download all activation, create 
            Emsr (str): name of EMSR activation in format EMSRXXX, ex. EMSR382
            Aoi (str): name of Area of interest in format AOIXX, ex. AOI02

        Returns:
            None
        """

        print("Download started...")
        rasterizer = RasterizeMask()

        # Use this only in case of DEBUG for future modification
        if utils.DEBUGDownload:
            print("DEBUG...")
            for key, value in self.dataDict.items():
                Emsr = key
                smoky_EMSR = ["EMSR369", "EMSR383", "EMSR387", "EMSR403", "EMSR408", "EMSR435", "EMSR541", "EMSR571", "EMSR578", "EMSR628", "EMSR632", "EMSR638", "EMSR524",  "EMSR534", "EMSR512", "EMSR545"]
                AOI_list = ["AOI01", "AOI02", "AOI03", "AOI04", "AOI05", "AOI06", "AOI07", "AOI08", "AOI09", "AOI10", "AOI11", "AOI12", "AOI13", "AOI14", "AOI15", "AOI16", "AOI17", "AOI18", "AOI19", "AOI20"]
                # AOI = "AOI01"
                start_date_ESMR = datetime.strptime(self.dataDict.get(Emsr), "%Y-%m-%dT%H:%M:%S.%f")
                end_date_ESMR = start_date_ESMR + timedelta(days = 30)

                if key in smoky_EMSR:
                    start_date_ESMR = start_date_ESMR + timedelta(days = 10)
                    end_date_ESMR = start_date_ESMR + timedelta(days = 50)

                for Aoi in AOI_list:
                    self.download_and_raster_images_EMSR(Emsr, 
                                                        Aoi, 
                                                        start_date_ESMR, 
                                                        end_date_ESMR,
                                                        pathFolder = utils.output_folder,
                                                        grading = grading, 
                                                        delineation = delineation, 
                                                        estimation = estimation,
                                                        )
            self.renameFiles()
            self.saveSatelliteDataCSV()

        # download all activation in output_folder
        elif allActivation:
            print("All activation downloading... 12 hour to complete")
            datasetPreconfigured = pd.read_csv(utils.datasetPreConfigured)
            for index, row in datasetPreconfigured.iterrows():
               
                row["folderPath"] = utils.output_folder + row["folderPath"]

                if not os.path.exists(row["folderPath"] + row["EMSR"] + "/" + row["AOI"] + "/"):

                    self.download_and_raster_images_EMSR(row["EMSR"], 
                                                        row["AOI"], 
                                                        row["suggested_startDate"], 
                                                        row["suggested_endDate"],
                                                        row["folderPath"],
                                                        grading = grading, 
                                                        delineation = delineation, 
                                                        estimation = estimation,
                                                        merged = True)
                    
                    

            print("Raster mask completed, writing satelliteDataCSV ... ")
            self.renameFiles()
            self.saveSatelliteDataCSV()

            print("Writing cloudCoverageCSV ...")
            open(utils.cloudCoveragecsv, 'w').close() #clean the file
            rasterizer.writeCloudcoverageCSV()
            
            print("Making landcovers...")
            # Land cover
            LC_download = DownloadLandCover()

            if "ESA" in utils.landCover:
                LC_download.download_ESA_WorldCover2020()
            if "LC9" in utils.landCover:
                LC_download.download_AnnualLandCover9()
            if "LC10" in utils.landCover:
                LC_download.download_Esri_LandCover10()

            print("Download Dataset done!!")

        # donwload a specific activation and area of interest
        elif Emsr != None and Aoi != None:

            print(f"Download {Emsr}_{Aoi}")
            datasetPreconfigured = pd.read_csv(utils.datasetPreConfigured)
            
            index = datasetPreconfigured.index[(datasetPreconfigured["EMSR"] == Emsr) & (datasetPreconfigured["AOI"] == Aoi) & (datasetPreconfigured["folderType"] != "subOptimal_smoke")].to_list()
            if bool(index) == False:
                print(f"{Emsr}_{Aoi} does not exists, wrong AOI or EMSR label")
                return
            datasetPreconfigured["folderPath"] = utils.output_folder + datasetPreconfigured["folderPath"].astype(str)
            path = datasetPreconfigured["folderPath"].loc[index].item() + datasetPreconfigured["EMSR"].loc[index].item() + "/" + datasetPreconfigured["AOI"].loc[index].item()
            
            # if the actiovation was previously downloaded, delete the folder and redownload it
            if os.path.exists(path):
                shutil.rmtree(path)

            self.download_and_raster_images_EMSR(   datasetPreconfigured["EMSR"].loc[index].item(), 
                                                    datasetPreconfigured["AOI"].loc[index].item(), 
                                                    datasetPreconfigured["suggested_startDate"].loc[index].item(), 
                                                    datasetPreconfigured["suggested_endDate"].loc[index].item(),
                                                    datasetPreconfigured["folderPath"].loc[index].item(),
                                                    grading = grading, 
                                                    delineation = delineation, 
                                                    estimation = estimation,
                                                    merged = True)

            print("Raster mask completed, writing satelliteDataCSV ... ")
            self.renameFiles(path)
            self.saveSatelliteDataCSV()

            open(utils.cloudCoveragecsv, 'w').close() #clean the file
            rasterizer.writeCloudcoverageCSV()
            
            print("Making landcovers...")
            # Land cover
            LC_download = DownloadLandCover()

            if "ESA" in utils.landCover:
                LC_download.download_ESA_WorldCover2020()
            if "LC9" in utils.landCover:
                LC_download.download_AnnualLandCover9()
            if "LC10" in utils.landCover:
                LC_download.download_Esri_LandCover10()
            
            print(f"Activation {Emsr}_{Aoi} done!!")
        
        else:
            print("Wrong configuration")

        return


    def retrieve_images_from_sentinel(self, path_output_folder: str, multipolygon_AOI: shapely.geometry.MultiPolygon, startDate: datetime, endDate: datetime, EMSR: str, AOI: str):
        
        """Download and retrieve the images of a given AOI with lowest cloud percentage in the given time interval

        Args:
            path_output_folder (str): output folder path
            multipolygon_AOI (shapely): multipolygon that defines boundaries of the AOI
            startDate (datetime): start date of the time interval of the query
            endDate (datetime): end date of the time interval of the query
            EMSR (str): activation name EMSR 
            AOI (str): relative area of interest

        Returns:
            files (List), acquistion (list): images and date of acquisition (today)
        """

        retr = ImageRetriever(utils.sh_client_id, utils.sh_client_secret)
        
        files, acquisition = retr.retrieve_images(evalscript = utils.copernicusFolderPath + utils.evalscript, 
                                                    data_collection = "",
                                                    geometry = multipolygon_AOI, 
                                                    min_dims = (utils.min_dims, utils.min_dims),
                                                    resolution = utils.resolution, 
                                                    start_date = startDate, 
                                                    end_date = endDate,
                                                    cache = True,
                                                    emsr = EMSR + "/" + AOI,
                                                    output_folder = path_output_folder) #utils.output_folder)
        
        return files, acquisition


    def download_and_raster_images_EMSR(self,
                                        EMSR: str, 
                                        AOI: str, 
                                        startDate: str, 
                                        endDate: str, 
                                        path_output_folder: str,
                                        grading:bool = True, 
                                        delineation:bool = True, 
                                        estimation: bool = True,
                                        merged: bool = False):
        
        """Download the images and create the relative raster mask for grading, delineation and first estimation 

        Args:
            EMSR (str): activation name EMSR 
            AOI (str): relative area of interest
            startDate (datetime): start date of the time interval of the query
            endDate (datetime): end date of the time interval of the query
            path_output_folder (str): output folder path
            grading (bool): create grading mask if relative json file exists
            delineation (bool): create delineation mask if relative json file exists 
            estimation (bool): create first estimation mask if relative json file exists 
            merged (bool): in case we want merge all images of the Aoi in one single image (NOTE: some Aoi take more than 20GB )
           
        Returns:
            None
        """

        
        if os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI):

            rasterizer = RasterizeMask()

            # DELINATION (DEL)                                    
            if delineation and os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/DEL/PRODUCT"):

                print(f"\nDownloading DELINATION - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read DEL_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/DEL/PRODUCT/" + EMSR + "_" + AOI + "_DEL_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]
                
                # create multipolygon of the Aoi
                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = MultiPolygon(polygon_AOI)

                if len(multipolygon_AOI) == 0:
                    print(f"{EMSR}/{AOI} Problem with Area of Interest Delineation -> possible rings as geometry? ")
                    file_object = open(utils.logPath, 'a')
                    file_object.write(EMSR + "/" + AOI + " Problem with Area of Interest Delineation -> possible rings as geometry? \n")
                    file_object.close()
                    return

                files, _ = self.retrieve_images_from_sentinel( path_output_folder, multipolygon_AOI, startDate, endDate, EMSR, AOI)

                if files is None:
                    return

                # create raster mask
                rasterizer.rasterizeMaskDelineation(output_folder = path_output_folder,
                                                    EMSR_AOI = EMSR + "/" + AOI + "/", 
                                                    files = files
                                                    )
                # merge in one single raster mask
                if merged:
                    rasterizer.rasterizeMaskDelineationMerged(  output_folder = path_output_folder,
                                                                EMSR_AOI = EMSR + "/" + AOI + "/", 
                                                                files = files
                                                                )

            # GRADING (GRA)
            if grading and os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/GRA/PRODUCT"):

                print(f"\nDownloading GRADING - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read GRA_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/GRA/PRODUCT/" + EMSR + "_" + AOI + "_GRA_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]

                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = MultiPolygon(polygon_AOI)

                if len(multipolygon_AOI) == 0:
                    print(f"{EMSR}/{AOI} Problem with Area of Interest Grading -> possible rings as geometry? ")
                    file_object = open(utils.logPath, 'a')
                    file_object.write(EMSR + "/" + AOI + " Problem with Area of Interest Grading -> possible rings as geometry? \n")
                    file_object.close()
                    return

                files, _ = self.retrieve_images_from_sentinel(path_output_folder, multipolygon_AOI, startDate, endDate, EMSR, AOI)
                
                if files is None:
                    return

                # create raster mask
                rasterizer.rasterizeMaskGrading(output_folder = path_output_folder,
                                                EMSR_AOI = EMSR + "/" + AOI + "/", 
                                                files = files
                                                )
                # merge in one single raster mask
                if merged:
                    rasterizer.rasterizeMaskGradingMerged(  output_folder = path_output_folder,
                                                            EMSR_AOI = EMSR + "/" + AOI + "/",
                                                            files = files 
                                                            )

            # FIRST ESTIMATION (FEP), only in case if there is no grading or delineation
            if estimation and os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/FEP/PRODUCT") and \
                not os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/GRA/PRODUCT") and \
                not os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/DEL/PRODUCT"):

                print(f"\nDownloading FIRST ESTIMATION - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read FEP_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/FEP/PRODUCT/" + EMSR + "_" + AOI + "_FEP_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]
                
                # NOTE: buffer(0) is a trick for fixing scenarios where polygons have overlapping coordinates 
                # create multipolygon of the Aoi
                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = MultiPolygon(polygon_AOI)

                if len(multipolygon_AOI) == 0:
                    print(EMSR + "/" + AOI + " Problem with Area of Interest Estimation -> possible rings as geometry? ")
                    file_object = open(utils.logPath, 'a')
                    file_object.write(EMSR + "/" + AOI + " Problem with Area of Interest Grading -> possible rings as geometry?  \n")
                    file_object.close()
                    return

                files, _ = self.retrieve_images_from_sentinel( path_output_folder, multipolygon_AOI, startDate, endDate, EMSR, AOI)
                
                if files is None:
                    return

                # create raster mask
                rasterizer.rasterizeMaskFEP(output_folder = path_output_folder,
                                            EMSR_AOI = EMSR + "/" + AOI + "/", 
                                            files = files
                                            )
                # merge in one single raster mask
                if merged:
                    rasterizer.rasterizeMaskFEPMerged(  output_folder = path_output_folder,
                                                        EMSR_AOI = EMSR + "/" + AOI + "/", 
                                                        files = files
                                                        )

            if not os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/GRA/PRODUCT") and \
                not os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/DEL/PRODUCT") and \
                not os.path.exists(utils.copernicusFolderPath + EMSR + "/" + AOI + "/FEP/PRODUCT"):
                print(f"=== {utils.copernicusFolderPath}{EMSR}/{AOI} === No DEL, GRA or FEP file found")
                file_object = open(utils.logPath, 'a')
                file_object.write(EMSR + "/" + AOI + " - No DEL, GRA or FEP file found \n")
                file_object.close()
                return
        else:
            print(f"Folder {utils.copernicusFolderPath + EMSR + '/' + AOI} does not exists")


    def saveSatelliteDataCSV(self):

        """Iterate over all dataset's folder and create satelliteData.csv with all information about each image, such as coordinate and resolution
        NOTE: due to low precision the resolution is read directly from tiff image meta_information and not from .json file"""
        
        # 15 decimals precision
        precision_decimals = 17
        getcontext().prec = precision_decimals    

        # clean the satelliteDataCSV file
        open(utils.satelliteData, 'w').close() 

        # Creation dataframe        
        columns_names = ["EMSR","AOI","folder","folderPath","activationDate","interval_startDate","interval_endDate","post_fire_acquisition","GRA","DEL","FEP","left_Long","bottom_Lat","right_Long","top_Lat","centerBoxLong","centerBoxLat","resolution_x","resolution_y","height","width","pixelBurned","country","koppen_group","koppen_subgroup"]
        satellite_data_df = pd.DataFrame(columns=columns_names)
        index = 0

        for dirpath, dirnames, filenames in os.walk(utils.output_folder):
            for filename in [f for f in filenames if "_S2L2A.tiff" in f]:

                if "dataColomba" in dirpath.split('/'):
                    optimal = "colomba"
                elif "dataSuboptimal" in dirpath.split('/'):
                    if "cloudy" in dirpath.split('/'):
                        optimal = "subOptimal_cloudy"
                    if "cloudyClean" in dirpath.split('/'):
                        optimal = "subOptimal_cloudyClean"
                    if "smoke" in dirpath.split('/'):
                        optimal = "subOptimal_smoke"
                    if "FEP" in dirpath.split('/'):
                        optimal = "subOptimal_FEP"                                 
                else:
                    optimal = "optimal"
                
                pathTiff = dirpath + "/" + filename
                pathJson = dirpath + "/" + filename.replace(".tiff", ".json")
                pathGRA = dirpath + "/" + filename.replace("_S2L2A.tiff", "_GRA.tif")
                pathDEL = dirpath + "/" + filename.replace("_S2L2A.tiff", "_DEL.tif")
                pathFEP = dirpath + "/" + filename.replace("_S2L2A.tiff", "_FEP.tif")

                # EMSR and AOI of image
                Emsr = filename.split("_")[0]
                Aoi = filename.split("_")[1]

                GRA_type = 0
                DEL_type = 0
                FEP_type = 0

                if os.path.exists(pathGRA):
                    GRA_type = 1
                    with rasterio.open(pathGRA) as gra_mask:
                        image_mask = gra_mask.read()
                        image_mask = np.squeeze(image_mask)
                        pixelBurned = np.count_nonzero(image_mask > 0)
                if os.path.exists(pathDEL):
                    DEL_type = 1
                    with rasterio.open(pathDEL) as del_mask:
                        image_mask = del_mask.read()
                        image_mask = np.squeeze(image_mask)
                        pixelBurned = np.count_nonzero(image_mask > 0)
                if os.path.exists(pathFEP):
                    FEP_type = 1
                    with rasterio.open(pathFEP) as fep_mask:
                        image_mask = fep_mask.read()
                        image_mask = np.squeeze(image_mask)
                        pixelBurned = np.count_nonzero(image_mask > 0)

                # resolution and bounding box
                bbox = np.zeros(4)
                with rasterio.open(pathTiff) as tif_image:
                    resolution_x = Decimal(tif_image.meta.get("transform")[0])
                    resolution_y = Decimal(tif_image.meta.get("transform")[4])
                    bbox[0] = Decimal(tif_image.profile["transform"][2])
                    bbox[1] = Decimal(tif_image.profile["transform"][5]) - (Decimal(tif_image.profile["height"]) * Decimal(abs(tif_image.profile["transform"][4])))
                    bbox[2] = Decimal(tif_image.profile["transform"][2]) + (Decimal(tif_image.profile["width"]) * Decimal(abs(tif_image.profile["transform"][0])))
                    bbox[3] = Decimal(tif_image.profile["transform"][5])
                
                # activation date of the event
                with open(utils.copernicusFolderPath + utils.dateJson) as json_date:
                    dataDict = json.load(json_date)
                    activationDate = dataDict.get(Emsr)

                # info from image json from sentinel hub
                with open(pathJson) as f:
                    json_file = json.load(f)
                    payload = json_file["payload"]
                    startDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["from"]
                    endDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["to"]
                    height = payload["output"]["height"]
                    width = payload["output"]["width"]

                    # acquisition date from sentinel2
                    if "dataSuboptimal" in dirpath.split('/'):
                        post_fire_acquisition = "Nan" # Not yet checked
                    else:
                        post_fire_acquisition = payload["acquisition_date"][0]

                
                # get country name where the event happened
                geolocator = Nominatim(user_agent="geoapiExercises")
                location = geolocator.reverse(f"{(Decimal(bbox[1])+Decimal(bbox[3]))/2}, {(Decimal(bbox[0])+Decimal(bbox[2]))/2}", language='en')
                address = location.raw['address']
                country = address.get('country', None)

                # get koppen climate classification label
                # Map downloaded from: https://figshare.com/articles/dataset/Present_and_future_K_ppen-Geiger_climate_classification_maps_at_1-km_resolution/6396959/2
                path_koppen_climate_map = "assets/koppen_climate/koppen_climate.tif"
                with rasterio.open(path_koppen_climate_map) as src:
                    row, col = src.index( float((Decimal(bbox[0])+Decimal(bbox[2]))/2), float((Decimal(bbox[1])+Decimal(bbox[3]))/2))
                    value_kop = src.read(1, window=((row, row+1), (col, col+1)))
                    value_kop = max(max(row) for row in value_kop)

                    if value_kop == 0:
                        value_kop = src.read(1, window=((row-5, row+5), (col-5, col+5)))
                        value_kop = max(max(row) for row in value_kop)

                        if value_kop == 0:
                            value_kop = src.read(1, window=((row-15, row+15), (col-15, col+15)))
                            value_kop = max(max(row) for row in value_kop)

                koppen_dict = { 0: "None", 1:  "Af", 2:  "Am", 3:  "Aw", 4:  "BWh", 5:  "BWk", 6:  "BSh", 7:  "BSk", 8:  "Csa", 9:  "Csb", 10: "Csc", 11: "Cwa", 12: "Cwb", 13: "Cwc", 14: "Cfa", 
                                15: "Cfb", 16: "Cfc", 17: "Dsa", 18: "Dsb", 19: "Dsc", 20: "Dsd", 21: "Dwa", 22: "Dwb", 23: "Dwc", 24: "Dwd", 25: "Dfa", 26: "Dfb", 27: "Dfc", 
                                28: "Dfd", 29: "ET", 30: "EF" }
                koppen_label = koppen_dict[value_kop]

                # store the data in dataframe
                satellite_data_df.loc[index] = [Emsr, Aoi, optimal, str(dirpath.replace(utils.output_folder, "")), str(activationDate), \
                                                str(startDate), str(endDate), str(post_fire_acquisition), GRA_type, DEL_type, FEP_type, \
                                                float(Decimal(bbox[0])), float(Decimal(bbox[1])), float(Decimal(bbox[2])), float(Decimal(bbox[3])), \
                                                float((Decimal(bbox[0])+Decimal(bbox[2]))/2), float((Decimal(bbox[1])+Decimal(bbox[3]))/2), float(resolution_x), float(resolution_y), \
                                                height, width, pixelBurned, country, koppen_label[0], koppen_label]
                index += 1

        # Sort and save the dataframe
        satellite_data_df_sorted = satellite_data_df.sort_values(['EMSR', 'AOI'], ascending = [True, True])
        satellite_data_df_sorted = satellite_data_df_sorted.round(decimals=precision_decimals)
        satellite_data_df_sorted.to_csv(utils.satelliteData, index=False)
            

    def downloadCloudCover(self):
        """Iterate over all dataset and download the cloud mask for each image"""

        for dirpath, dirnames, filenames in os.walk(utils.output_folder):
            for image in [f for f in filenames if "_S2L2A.tiff" in f]:
                
                print(dirpath)
                rasterizer = RasterizeMask()
                _, _ = rasterizer.createCloudMaskfromFilePath(image=image, path=dirpath)
    

    def downloadLandCover(self, redownload: bool = False):
        """Iterate over all dataset and download the land cover for each image"""
        
        if redownload == True:
            for dirpath, _, filenames in os.walk(utils.output_folder):
                for image in [f for f in filenames if "_LC.tif" in f]:
                    path = dirpath + "/" + image
                    pathRGB = path.replace("_LC.tif", "_LC.png")
                    if os.path.exists(path):
                        os.remove(path)
                    if os.path.exists(pathRGB):
                        os.remove(pathRGB)

        LC_download = DownloadLandCover()

        if "ESA" in utils.landCover:
            LC_download.download_ESA_WorldCover2020()
        if "LC9" in utils.landCover:
            LC_download.download_AnnualLandCover9()
        if "LC10" in utils.landCover:
            LC_download.download_Esri_LandCover10()


    def renameFiles(self, pathToSingleEMSR: str = ""):
        """Rename hash code with respective EMSR and AOI number"""

        # All activations
        if pathToSingleEMSR == "":
            for dirpath, _, filenames in os.walk(utils.output_folder):
                for filename in [f for f in filenames if "_S2L2A.png" in f]:

                    renamedTiff = filename.split(".")[0] + ".tiff"
                    renamedJson = filename.split(".")[0] + ".json"
                    renamedFolder = dirpath.replace(dirpath.split("/")[-1], filename.split(".")[0].replace("_S2L2A", ""))
                    
                    if os.path.exists(dirpath + "/response.tiff"):
                        os.rename(dirpath + "/response.tiff", dirpath + "/" + renamedTiff)
                        os.rename(dirpath + "/request.json", dirpath + "/" + renamedJson)
                        os.rename(dirpath, renamedFolder)

        # Single activation
        else:
            for dirpath, _, filenames in os.walk(pathToSingleEMSR):
                for filename in [f for f in filenames if "_S2L2A.png" in f]:

                    renamedTiff = filename.split(".")[0] + ".tiff"
                    renamedJson = filename.split(".")[0] + ".json"
                    renamedFolder = dirpath.replace(dirpath.split("/")[-1], filename.split(".")[0].replace("_S2L2A", ""))
                    
                    if os.path.exists(dirpath + "/response.tiff"):
                        os.rename(dirpath + "/response.tiff", dirpath + "/" + renamedTiff)
                        os.rename(dirpath + "/request.json", dirpath + "/" + renamedJson)
                        os.rename(dirpath, renamedFolder)
    

    def copyAllImageToRGBFolder(self):
        """Iterate over all dataset and copy each image .png to RGB_folder"""

        for root, dirs, files in os.walk(utils.output_folder):
            for file in files:
                if file.endswith(".png"):
                    
                    path_file = os.path.join(root,file)
                    path_dest = root.replace(utils.output_folder,utils.RGB_folder)

                    if not os.path.isdir(path_dest):
                        os.makedirs(path_dest)

                    shutil.copy2(path_file, path_dest) 
    

    def createPreconfiguredCSV(self):
        """Iterate over satelliteData.csv to retrieve all the date information about the images location in the dataset. 
        Then save all the information for future replication of the dataset. 
        NOTE: USE THIS FUNCTION ONLY IN CASE OF DEBUG
        """
 
        open(utils.datasetPreConfigured, 'w').close() #clean the file
        file_object = open(utils.datasetPreConfigured, 'a')
        file_object.write("EMSR,AOI,folderType,folderPath,activationDate,suggested_startDate,suggested_endDate\n")
        file_object.close()

        df = pd.read_csv(utils.satelliteData)
        df = df.drop(columns=["GRA","DEL","FEP", "left_Long", "bottom_Lat", "right_Long", "top_Lat", "centerBoxLong", "centerBoxLat", "height", "width", "resolution_x", "resolution_y"])

        df["folderPath"] = df["folderPath"].apply(lambda pathFolder: "/".join(pathFolder.split("/")[:2]) if pathFolder.split("/")[0] == "dataSuboptimal" else pathFolder.split("/")[0])
        df["folderPath"] = df["folderPath"].apply(lambda pathFolder: pathFolder+"/")
        # df = df.groupby(["EMSR", "AOI"], as_index = False)
        df = df.drop_duplicates()

        df = df.sort_values(['EMSR', 'AOI'], ascending = [True, True])
        df.to_csv(utils.datasetPreConfigured, mode='a', index=False, header=False)
