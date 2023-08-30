from datetime import datetime, timedelta
from typing import List, Tuple, Union
from pathlib import Path
import datetime
import json

import numpy as np
import sentinelhub as sh
from geojson_pydantic import MultiPolygon, Polygon
from shapely.geometry import shape

import src.utils_variable as utils

import warnings
warnings.filterwarnings('ignore')


class ImageRetriever():

    def __init__(self, id: str, secret: str) -> None:
        self.sh_config = sh.SHConfig(hide_credentials=True)
        self.sh_config.sh_client_id = id
        self.sh_config.sh_client_secret = secret
        self.catalog = sh.SentinelHubCatalog(config=self.sh_config)
        # self.date_format = "%Y-%m-%d"

    
    def _check_geometry_area(self, bboxes: List[sh.BBox], resolution: int, emsr: str):
        """verify if the input geometry area is lower than the maximum allowed

        Args:
            bboxes (List[sh.BBox]): input bboxes
            resolution (int): resolution

        Raises:
            ValidationError: validation error
        """
        area = 0
        for bbox in bboxes:
            w, h = sh.bbox_to_dimensions(bbox, resolution=resolution)
            area += w * h
        max_w, max_h = (15000,15000)
        max_area = max_w * max_h
        if area > max_area:
            print("Area exceeded")
            file_object = open(utils.logPath, 'a')
            file_object.write(emsr + " - Area exceeded\n")
            file_object.close()
            pass


    def retrieve_images(self, evalscript: str, data_collection: str, geometry: Union[Polygon, MultiPolygon],
                        min_dims: Tuple[int, int], resolution: int, start_date: datetime, end_date: datetime,
                        cache: bool, emsr: str, output_folder: str) -> list:
        """obtain images from SentinelHub

        Args:
            evalscript (str): evalscript file name
            data_collection (str): data collection string according to sh.DataCollection
            geometry (Union[Polygon, MultiPolygon]): geometry
            min_dims (Tuple[int, int]): minimum dimensions of the AoI
            resolution (int): resolution
            start_date (datetime): start date
            end_date (datetime): end date
            cache(bool): is caching is used

        Returns:
            list: list of results
        """

        bboxes = self.get_area_splits(geometry, min_dims, resolution)

        # data_collection = sh.DataCollection[data_collection]
        data_collection=sh.DataCollection.SENTINEL2_L2A
        mosaicking = None
        query = None
        files = []
        acquisition_dates = []

        # INFO: https://docs.sentinel-hub.com/api/latest/data/sentinel-2-l2a/
        if data_collection in (sh.DataCollection.SENTINEL2_L1C, sh.DataCollection.SENTINEL2_L2A):
            mosaicking = "leastCC" # leastRecent leastCC
            query = {"eo:cloud_cover": {"lt": int(utils.cloudCoverMax * 100)}}
        # evalscript = (self.config.evalscript_dir / evalscript).read_text()
        evalscript =  Path(evalscript).read_text() #Path("evalScriptV3.js").read_text()
        try:
            self._check_geometry_area(bboxes, resolution, emsr)
            for box in bboxes:
                file, dates = self.download_images(evalscript=evalscript,
                                                   data_collection=data_collection,
                                                   resolution=resolution,
                                                   start_date=start_date,
                                                   end_date=end_date,
                                                   bbox=box,
                                                   maxcc=utils.cloudCoverMax,
                                                   mosaicking_order=mosaicking,
                                                   query=query,
                                                   cache=cache,
                                                   emsr=emsr,
                                                   output_folder=output_folder)
                if file is None:
                    return None, None
                    
                files.append(file)
                acquisition_dates.append([d.utcnow().isoformat() for d in dates])

            return files, acquisition_dates

        except Exception as te:
            print(str(te))
            pass

  
    def download_images(self,
                        evalscript: str,
                        data_collection: sh.DataCollection,
                        resolution: int,
                        start_date: str,
                        end_date: str,
                        bbox: sh.BBox,
                        maxcc: float,
                        cache: bool,
                        mosaicking_order: str = None,
                        query: dict = None,
                        emsr: str = None,
                        output_folder: str = "data/"):
        """Perform the request to Sentinel Hub APIs

        Args:
            evalscript (str): script containing the collection configuration
            data_collection (sh.DataCollection): Sentinel data collection
            resolution (int): resolution
            start_date (str): start date
            end_date (str): end date
            bbox (sh.BBox): bbox
            maxcc(float): maxcc allowed
            cache(bool): if caching is used
            query(dict): query for search activations
            mosaicking_order(str): mosaicking order (for Sentinel2 only)
        """

        time_difference = timedelta(days=1)
        time_interval = start_date, end_date

        search_iterator = self.catalog.search(data_collection, bbox=bbox, time=time_interval, query=query)
        all_timestamps = search_iterator.get_timestamps()
        unique_acquisitions = sh.filter_times(all_timestamps, time_difference)

        # print(all_timestamps)

        if len(unique_acquisitions) == 0:
            print(f"No images found in period {start_date} - {end_date}.")
            file_object = open(utils.logPath, 'a')
            file_object.write(f"{emsr} No images found in time interval {start_date} - {end_date} \n")
            file_object.close()
            return None, None

        # self.log.debug("Downloading %s, from %s to %s, area %s", data_collection, start_date, end_date, bbox)
        request = sh.SentinelHubRequest(evalscript = evalscript,
                                        data_folder = output_folder+emsr,
                                        input_data=[
                                            sh.SentinelHubRequest.input_data(data_collection=data_collection,
                                                                             time_interval=(start_date, end_date),
                                                                             maxcc=maxcc,
                                                                             mosaicking_order=mosaicking_order)
                                        ],
                                        responses=[sh.SentinelHubRequest.output_response("default", sh.MimeType.TIFF)],
                                        bbox=bbox,
                                        size=sh.bbox_to_dimensions(bbox, resolution=resolution),
                                        config=self.sh_config
                                        )
        # with this kind of process, just one file per request is returned
        # raise_download_errors will trigger DownloadFailedException on failure
        try:
            request.get_data(save_data=True, max_threads=4, raise_download_errors=True, redownload=(not cache))
        except Exception as e:
            print(f"Cannot download images: {str(e)}")
            pass

        file_name = request.get_filename_list()[0]

        # save the acquisition date inside json
        path_to_json = Path(output_folder+emsr+"/"+file_name.replace("response.tiff","request.json"))

        with open(path_to_json, "r") as f:
            json_request = json.load(f)
        
        payload = json_request['payload']
        payload['acquisition_date'] = [acquisition.strftime('%Y/%m/%d_%H:%M:%S') for acquisition in unique_acquisitions]
        json_request['payload'] = payload

        with open(path_to_json, 'w') as f:
            json.dump(json_request, f)

        return file_name, unique_acquisitions

  
    def get_area_splits(self, geometry: Union[Polygon, MultiPolygon], min_dims: Tuple[int, int], resolution: int):
        """        Retrieves the bounding box of the region of interest starting from polygon coordinates. If the area is too large,
        it is split into smaller areas. Returns the list of the bounding boxes and a tuple containing the widths and heights
        (in pixels) of the images of the areas, together with the specified resolution
        Args:
            geometry (Union[Polygon, MultiPolygon]): input polygon or multipolygon, only the bounds will be used.
            min_dims (Tuple[int, int]): minimum dimensions, in pixels, of the bounding box
            resolution (int): resolution
        Returns:
            List[sh.BBox]: list of bounding boxes derived from the splits
        """

        bbox = sh.BBox(bbox=shape(geometry).bounds, crs=sh.CRS.WGS84)
        width, height = sh.bbox_to_dimensions(bbox, resolution=resolution)
        # ugly workaround: whenever the box dimensions are lower than the minimum tiling size,
        # buffer the original geometry until compliant.
        h, w = min_dims
        while width < w or height < h:
            geometry = shape(geometry).buffer(0.01)
            width, height = sh.bbox_to_dimensions(sh.BBox(bbox=shape(geometry).bounds, crs=sh.CRS.WGS84),
                                                  resolution=resolution)
        # once the minimum requirements are satistified, produce a series of tiles
        # computing the number of vertical and horizontal splits, based on the _max_ dimensions
        h_splits = int(np.ceil(width / utils.splitMaxDimension))
        v_splits = int(np.ceil(height / utils.splitMaxDimension))

        bbox_splitter = sh.BBoxSplitter([shape(geometry)], sh.CRS.WGS84, (h_splits, v_splits))
        return bbox_splitter.get_bbox_list()