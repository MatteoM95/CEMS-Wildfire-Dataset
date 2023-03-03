# CEMS-Wildfire-Dataset
Copernicus emergency management service wildfire dataset from 2019 to 2022.

# Introduction
This dataset wiil be part of European project [SAFERS](https://safers-project.eu/) ‘Structured Approaches for Forest Fire Emergencies in Resilient Societies’. It is going to create an open and integrated platform featuring a forest fire Decision Support System. The platform will use information from different sources: earth observations from Copernicus and GEOSS, fire sensors in forests, topographic data, weather forecasts and even crowdsourced data from social media and other apps that can be used by citizens and first responders to provide situational in-field information.
A sample from dataset is available [here](ZZZZ)


## Installation packages
First install all requirements in [assets/requirements.txt](ZZZZZ) and cloudSen12 library
```
pip install -r /assets/requirements.txt
pip install cloudsen12
```
In case error with cartopy library, install
```
sudo apt -y install libgeos-dev
```


## How to run and create the dataset
```python
"""NOTE: All variable and parameter for customization such as path/to/folder are in src.utils_variables"""
from src.downloadDataset import DownloadCEMSDataset as CemsDataset

import warnings
warnings.filterwarnings('ignore')

download_cems = CemsDataset()

"""Uncomment this line if you would download and create the whole dataset"""
# download_cems.download_EMSR_Manager(allActivation = True)

"""Uncomment this line if you would download a specif activation, pay attention to the format: EMSRXXX and AOIXX"""
# download_cems.download_EMSR_Manager(Emsr = "EMSR638", Aoi = "AOI01")

"""Uncomment this line if you would create/recreate the cloud cover for the whole dataset"""
# download_cems.downloadCloudCover()

"""Uncomment this line if you would re/download land cover for the whole dataset"""
# download_cems.downloadLandCover()

"""Uncomment this line if you would re/save satellite data with all information about the dataset"""
# download_cems.saveSatelliteDataCSV()

"""Uncomment this line if you have changed the folder arrangement, so it will be possible replicate the dataset arrangement
NOTE: USE THIS FUNCTION ONLY IN CASE OF DEBUG"""
# download_cems.createPreconfiguredCSV()

"""Uncomment this line if you would copy all .png images in a separate folder"""
# download_cems.copyAllImageToRGBFolder()
```


## Landcovers

In addition to wildfire delineation mask and clod masks, also landcovers for each image are create. In particular the model considered are: [**ESRI 10m Annual Land Use Land Cover (2017-2021)**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1), [**ESRI 2020 Global Land Use Land Cover**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1), [**ESA WorldCover 10 m 2020**](https://esa-worldcover.org/en/data-access).
All landcover models are based on sentinel2 10-meter resolution images.

#### ESRI 10m Annual Land Use Land Cover (2017-2021) (9 classes)
There are 9 classes:

Label number | Class | Class definitions | Color
--- | --- | ---
1 | `Water` | Water Areas where water was predominantly present throughout the year; may not cover areas with sporadic or ephemeral water; contains little to no sparse vegetation, no rock outcrop nor built up features like docks; examples: rivers, ponds, lakes, oceans, flooded salt plains. | ![#1A5BAB](https://placehold.co/15x15/1A5BAB/1A5BAB.png)
2 | `Trees` | Any significant clustering of tall (~15 feet or higher) dense vegetation, typically with a closed or dense canopy; examples: wooded vegetation, clusters of dense tall vegetation within savannas, plantations, swamp or mangroves (dense/tall vegetation with ephemeral water or canopy too thick to detect water underneath). | ![#358221](https://placehold.co/15x15/358221/358221.png)
4 | `Flooded Vegetation` | Areas of any type of vegetation with obvious intermixing of water throughout a majority of the year; seasonally flooded area that is a mix of grass/shrub/trees/bare ground; examples: flooded mangroves, emergent vegetation, rice paddies and other heavily irrigated and inundated agriculture. | ![#87D19E](https://placehold.co/15x15/87D19E/87D19E.png)
5 | `Crops` | Human planted/plotted cereals, grasses, and crops not at tree height; examples: corn, wheat, soy, fallow plots of structured land. | ![#FFDB5C](https://placehold.co/15x15/FFDB5C/FFDB5C.png)
7 | `Built Area` |Human made structures; major road and rail networks; large homogeneous impervious surfaces including parking structures, office buildings and residential housing; examples: houses, dense villages / towns / cities, paved roads, asphalt. | ![#ED022A](https://placehold.co/15x15/ED022A/ED022A.png)
8 | `Bare Ground` | Areas of rock or soil with very sparse to no vegetation for the entire year; large areas of sand and deserts with no to little vegetation; examples: exposed rock or soil, desert and sand dunes, dry salt flats/pans, dried lake beds, mines. | ![#EDE9E4](https://placehold.co/15x15/EDE9E4/EDE9E4.png)
9 | `Snow/Ice` | Large homogeneous areas of permanent snow or ice, typically only in mountain areas or highest latitudes; examples: glaciers, permanent snowpack, snow fields. | ![#F2FAFF](https://placehold.co/15x15/F2FAFF/F2FAFF.png)
10 | `Clouds` | No land cover information due to persistent cloud cover. | ![#C8C8C8](https://placehold.co/15x15/C8C8C8/C8C8C8.png)
11 | `Rangeland` | Open areas covered in homogeneous grasses with little to no taller vegetation; wild cereals and grasses with no obvious human plotting (i.e., not a plotted field); examples: natural meadows and fields with sparse to no tree cover, open savanna with few to no trees, parks/golf courses/lawns, pastures. | ![#C6AD8D](https://placehold.co/15x15/C6AD8D/C6AD8D.png)

#### ESRI 2020 Global Land Use Land Cover


#### ESA WorldCover 10m 2020

[Esa worldcover manual](https://esa-worldcover.s3.amazonaws.com/v100/2020/docs/WorldCover_PUM_V1.0.pdf)

in case broken git: git reset --hard origin/master
