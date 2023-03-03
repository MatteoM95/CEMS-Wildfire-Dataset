import pandas as pd
import numpy as np
import json
import os
import shutil

from datetime import datetime, timedelta
from decimal import *
import rasterio.features

import shapely
from shapely.geometry import shape, GeometryCollection

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


    def download_EMSR_Manager(self, allActivation: bool = False, Emsr: str = None, Aoi: str = None):
        
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

        # Use this only in case of download activation and make future modification
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
                                                        pathFolder = utils.output_folder
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
                                                        merged = True)
            print("Raster mask completed")
            self.renameFiles()
            self.saveSatelliteDataCSV()

            open(utils.cloudCoveragecsv, 'w').close() #clean the file
            rasterizer.writeCloudcoverageCSV()

            # Land cover
            LC_download = DownloadLandCover()

            if "ESA" in utils.landCover:
                LC_download.download_ESA_WorldCover2020()
            if "LC9" in utils.landCover:
                LC_download.download_AnnualLandCover9()
            if "LC10" in utils.landCover:
                LC_download.download_Esri_LandCover10()

            print("Download Dataset done!!")

        # donwload a specif activation and area of interest
        elif Emsr != None and Aoi != None:

            print(f"Download {Emsr}_{Aoi}")
            datasetPreconfigured = pd.read_csv(utils.datasetPreConfigured)
            
            index = datasetPreconfigured.index[(datasetPreconfigured["EMSR"] == Emsr) & (datasetPreconfigured["AOI"] == Aoi) & (datasetPreconfigured["folderType"] != "subOptimal_smoke")].to_list()
            datasetPreconfigured["folderPath"] = utils.output_folder + datasetPreconfigured["folderPath"].astype(str)
            path = datasetPreconfigured["folderPath"].loc[index].item() + datasetPreconfigured["EMSR"].loc[index].item() + "/" + datasetPreconfigured["AOI"].loc[index].item()
            
            # if the actiovation was already downloaded, delete the folder and redownload
            if os.path.exists(path):
                shutil.rmtree(path)
            self.download_and_raster_images_EMSR(datasetPreconfigured["EMSR"].loc[index].item(), 
                                                datasetPreconfigured["AOI"].loc[index].item(), 
                                                datasetPreconfigured["suggested_startDate"].loc[index].item(), 
                                                datasetPreconfigured["suggested_endDate"].loc[index].item(),
                                                datasetPreconfigured["folderPath"].loc[index].item(),
                                                merged = True)

            print("Raster mask completed")
            self.renameFiles(path)
            self.saveSatelliteDataCSV()

            open(utils.cloudCoveragecsv, 'w').close() #clean the file
            rasterizer.writeCloudcoverageCSV()

            # Land cover
            LC_download = DownloadLandCover()

            if  "ESA" in utils.landCover:
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

                print(f"Downloading DELINATION - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read DEL_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/DEL/PRODUCT/" + EMSR + "_" + AOI + "_DEL_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]
                
                # create multipolygon of the Aoi
                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = shapely.geometry.MultiPolygon(polygon_AOI)

                if len(multipolygon_AOI) == 0:
                    print(f"{EMSR}/{AOI} Problem with Area of Interest Grading -> possible rings as geometry? ")
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

                print(f"Downloading GRADING - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read GRA_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/GRA/PRODUCT/" + EMSR + "_" + AOI + "_GRA_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]
                
                # NOTE: buffer(0) is a trick for fixing scenarios where polygons have overlapping coordinates 
                # create multipolygon of the Aoi
                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = shapely.geometry.MultiPolygon(polygon_AOI)

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

                print(f"Downloading FIRST ESTIMATION - {EMSR}_{AOI} - Path output:  {path_output_folder} - time: {datetime.now()} " )

                # read FEP_PRODUCT_areaOfInterestA.json
                with open(utils.copernicusFolderPath + EMSR + "/" + AOI + "/FEP/PRODUCT/" + EMSR + "_" + AOI + "_FEP_PRODUCT_areaOfInterestA.json") as f:
                    features = json.load(f)["features"]
                
                # NOTE: buffer(0) is a trick for fixing scenarios where polygons have overlapping coordinates 
                # create multipolygon of the Aoi
                polygon_AOI = GeometryCollection([shape(feature["geometry"]).buffer(0) for feature in features if "type" in list(feature["geometry"]) and feature["geometry"]["type"] == "Polygon"  ]  )
                multipolygon_AOI = shapely.geometry.MultiPolygon(polygon_AOI)

                if len(multipolygon_AOI) == 0:
                    print(EMSR + "/" + AOI + " Problem with Area of Interest Grading -> possible rings as geometry? ")
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


    def saveSatelliteDataCSV(self):

        """Iterate over all dataset's folder and create satelliteData.csv with all information about each image, such as coordinate and resolution
        NOTE: due to low precision the resolution is read directly from tiff image meta_information and not from .json file"""
        
        # 20 decimals precision
        getcontext().prec = 20      

        #clean the file
        open(utils.satelliteData, 'w').close() 

        #Creation CSV
        if os.stat(utils.satelliteData).st_size == 0:
            file_object = open(utils.satelliteData, 'a')
            file_object.write("EMSR,AOI,folder,folderPath,activationDate,interval_startDate,interval_endDate,GRA,DEL,FEP,left_Long,bottom_Lat,right_Long,top_Lat,centerBoxLong,centerBoxLat,height,width,resolution_x,resolution_y\n")
            file_object.close()

        for dirpath, dirnames, filenames in os.walk(utils.output_folder):
            for filename in [f for f in filenames if "_S2L2A.tiff" in f]:

                if "dataSuboptimal" in dirpath.split('/'):
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

                Emsr = filename.split("_")[0]
                Aoi = filename.split("_")[1]

                GRA_type = 0
                DEL_type = 0
                FEP_type = 0

                if os.path.exists(pathGRA):
                    GRA_type = 1
                if os.path.exists(pathDEL):
                    DEL_type = 1
                if os.path.exists(pathFEP):
                    FEP_type = 1

                # resolution and bounding box
                bbox = np.zeros(4)
                with rasterio.open(pathTiff) as tif_image:
                    resolution_x = Decimal(tif_image.meta.get("transform")[0])
                    resolution_y = Decimal(tif_image.meta.get("transform")[4])
                    bbox[0] = Decimal(tif_image.profile["transform"][2])
                    bbox[1] = Decimal(tif_image.profile["transform"][5]) - (Decimal(tif_image.profile["height"]) * Decimal(abs(tif_image.profile["transform"][4])))
                    bbox[2] = Decimal(tif_image.profile["transform"][2]) + (Decimal(tif_image.profile["width"]) * Decimal(abs(tif_image.profile["transform"][0])))
                    bbox[3] = Decimal(tif_image.profile["transform"][5])

                # activation date
                with open(utils.copernicusFolderPath + utils.dateJson) as json_date:
                    dataDict = json.load(json_date)
                    activationDate = dataDict.get(Emsr)

                # from json
                with open(pathJson) as f:
                    json_file = json.load(f)
                    payload = json_file["payload"]
                    startDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["from"]
                    endDate = payload["input"]["data"][0]["dataFilter"]["timeRange"]["to"]
                    height = payload["output"]["height"]
                    width = payload["output"]["width"]
                
                # line string
                CSV_line = Emsr + "," + Aoi + "," + optimal + "," + str(dirpath.replace(utils.output_folder, "")) + "," + str(activationDate) + "," \
                        + str(startDate) + "," + str(endDate) + "," + str(GRA_type) + "," + str(DEL_type) + "," + str(FEP_type) + "," \
                        + str(Decimal(bbox[0])) + "," + str(Decimal(bbox[1])) + "," + str(Decimal(bbox[2])) + "," + str(Decimal(bbox[3])) + "," \
                        + str((Decimal(bbox[0])+Decimal(bbox[2]))/2) + "," + str((Decimal(bbox[1])+Decimal(bbox[3]))/2) + "," + str(height) + "," \
                        + str(width) + "," + str(resolution_x) + "," + str(resolution_y) + "\n"

                file_object = open(utils.satelliteData, 'a')
                file_object.write(CSV_line)
                file_object.close()

        # Sort and order the data
        df = pd.read_csv(utils.satelliteData, float_precision='round_trip')

        # clean the file
        open(utils.satelliteData, 'w').close()

        # write headline
        file_object = open(utils.satelliteData, 'a')
        file_object.write("EMSR,AOI,folder,folderPath,activationDate,interval_startDate,interval_endDate,GRA,DEL,FEP,left_Long,bottom_Lat,right_Long,top_Lat,centerBoxLong,centerBoxLat,height,width,resolution_x,resolution_y\n")
        file_object.close()

        # save sorted dataframe
        df = df.sort_values(['EMSR', 'AOI'], ascending = [True, True])
        df.to_csv(utils.satelliteData, mode='a', index=False, header=False, float_format='%.20f')


    def downloadCloudCover(self):
        """Iterate over all dataset and download the cloud mask for each image"""

        for dirpath, dirnames, filenames in os.walk(utils.output_folder):
            for image in [f for f in filenames if "_S2L2A.tiff" in f]:
                
                print(dirpath)
                rasterizer = RasterizeMask()
                _, _ = rasterizer.createCloudMaskfromFilePath(image=image, path=dirpath)
    

    def downloadLandCover(self):
        """Iterate over all dataset and download the land cover for each image"""
        
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
        




















        
# To find Sentinel-2 images that are similar to a reference image over a certain time interval, 
# you could use a combination of the approaches I mentioned earlier:

# 1 Use a spectral matching algorithm: This approach involves comparing the spectral characteristics of the reference 
# image with those of other Sentinel-2 images. You could use a spectral matching algorithm to identify images with similar 
# spectral profiles to the reference image.

# 2 Use image registration: Image registration is the process of aligning two or more images of the same scene so 
# that they can be compared. You could use image registration to align the reference image with other Sentinel-2 images
# and then compare them to identify those that are most similar to the reference image.

# 3 Use image segmentation: Image segmentation is the process of dividing an image into distinct regions or segments based 
# on certain characteristics. You could use image segmentation to divide the reference image into segments and then compare 
# the segments with those in other Sentinel-2 images to identify those that are most similar.

# 4 Use a machine learning model: You could train a machine learning model to identify Sentinel-2 images that are 
# similar to the reference image. To do this, you would need to have a large dataset of Sentinel-2 images along with 
# labels indicating which images are similar to the reference image. The model could then be used to predict which of the 
# other Sentinel-2 images are most similar to the reference image.

# You could also try using a combination of these approaches to improve the accuracy of your results.
# To find the Sentinel-2 images that you want to download, you can use the Sentinel Hub EO Browser or the Sentinel-2 Toolbox. 
# These tools allow you to search for Sentinel-2 images based on various criteria, including the time interval you are interested in. 
# You can then use the tools to download the images that you have identified as being most similar to your reference image.
    