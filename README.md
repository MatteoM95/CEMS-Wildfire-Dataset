# Copernicus Emergency Management Service - Wildfire Dataset (2017- 2023)
The Copernicus Emergency Management Service (CEMS) Wildfire dataset spans from June 2017 to April 2023. The dataset includes Sentinel-2 images related to wildfires, along with their respective severity and delineation masks. Additionally, the dataset is enhanced with cloud and landcover masks, providing more valuable information for future training of a semantic segmentation model. The dataset comprises over 500+ high-quality images, suitable for subsequent semantic segmentation model training. The dataset is available on [**Huggingface**](https://huggingface.co/datasets/links-ads/wildfires-cems)


## Introduction and Data sources

The **Copernicus Emergency Management Service (CEMS)** is a component of the European Union's Copernicus Programme. It provides rapid geospatial information during emergencies, and damage assessment for events such as floods, wildfires, and earthquakes. In particular, [**Copernicus Rapid Mapping**](https://emergency.copernicus.eu/mapping/list-of-activations-rapid) provides on-demand mapping services in cases of various natural disasters, offering detailed and up-to-date geospatial information that assists in disaster management and risk assessment. 

The satellite imagery comes from Sentinel-2 L2A, which spans across 12 distinct bands of the light spectrum, each with a resolution of 10 meters. The Sentinel Level-2A data undergoes an atmospheric correction process to adjust for the reflectance values influenced by the moisture in the atmosphere. These images are downloaded using SentinelHub APIs. 

## Data structure and organization

The structure of the dataset is the following:

```
dataset/
├── dataOptimal/
│   ├── EMSRXXX/
│   │   ├── AOIYY/
│   │   │   ├── EMSRXXX_AOIYY_01/
│   │   │   │   ├── EMSRXXX_AOIYY_01_Annual9_LC.png           # Landcover data Lulc (9 classes)
│   │   │   │   ├── EMSRXXX_AOIYY_01_Annual9_LC.tif           # Landcover data Lulc (9 classes) in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_CM.png                   # Cloud mask generated from cloudSen12 
│   │   │   │   ├── EMSRXXX_AOIYY_01_CM.tif                   # Cloud mask from cloudSen12 in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_DEL.png                  # Delineation mask
│   │   │   │   ├── EMSRXXX_AOIYY_01_DEL.tif                  # Delineation mask in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_ESA_LC.png               # Landcover data ESA WorldCover 2020
│   │   │   │   ├── EMSRXXX_AOIYY_01_ESA_LC.tif               # Landcover data ESA WorldCover 2020 in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_Esri10_LC.png            # Landcover data 2020 Global 10 Class (LULC)
│   │   │   │   ├── EMSRXXX_AOIYY_01_Esri10_LC.tif            # Landcover data 2020 Global 10 Class (LULC) in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_GRA.png                  # Grading mask 
│   │   │   │   ├── EMSRXXX_AOIYY_01_GRA.tif                  # Grading mask in georeferenced format
│   │   │   │   ├── EMSRXXX_AOIYY_01_S2L2A.json               # Image additional metadata from SentinelHub
│   │   │   │   ├── EMSRXXX_AOIYY_01_S2L2A.png                # Sentinel2 image
│   │   │   │   └── EMSRXXX_AOIYY_01_S2L2A.tiff               # Sentinel2 image in georeferenced format
│   │   │   └── EMSRXXX_AOIYY_01_merged.png                   # Merge between several tiles from sentinelHub
│   │   │  
│   │   ├── EMSRXXX_AOIYY_02/
│   │   │   └── ...
│   │
├── dataSuboptimal/
│   └── ...

```

A [sample](https://github.com/MatteoM95/CEMS-Wildfire-Dataset/tree/main/assets/sample/EMSR382/AOI01) from the dataset is made available to give you a representative overview of the data structure and accompanying metadata.


All the informations of the dataset are available inside [csv_files/](https://github.com/MatteoM95/CEMS-Wildfire-Dataset/tree/main/csv_files) folder:
- *dataset_Preconfigured.csv*: It contains all the activations from the dataset, including the activation date of the event and the interval date for SentinelHub API
- *satelliteData.csv*: All information about each image is stored here. 
- *log.txt*: general log for errors and messages.






## Installation packages
First, install all requirements in [assets/requirements.txt](https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/requirements.txt) and cloudSen12 library
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
# download_cems.download_EMSR_Manager(Emsr = "EMSR382", Aoi = "AOI01", grading = True, delineation = True, estimation = True)

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




## Wildfire masks

From Copernicus Rapid Mapping for each activation under the tag **Wildfire** are available different post-fire products:
- **FEP** (First Estimation): a first estimation of the burned area.
- **DEL** (Delineation): a delineation of the area affected by the wildfire.
- **GRA** (Grading): a detailed description about the severity of the burned area.

Each product includes metadata and associated JSON files which contain geographical details about the affected areas.

NOTE: Since I do not have the license to distribute those files then they must be retrived directly from Copenicus website. I left only the folder structure in copernicusData folder.

On Copernicus site are available several georeferenced data in GeoJSON format. For this dataset are considered only those files for each activation that are formatted with the following string:
- ``EMSRXXX_AOIYY_TYPE_PRODUCT_areaOfInterestA.json``, where TYPE can be `GRA`, `DEL` or `FEP`, defines the AOI YY of that particular activation
XXX where the event happened.
- ``EMSRXXX_AOIYY_TYPE_PRODUCT_observedEventA.json``, where TYPE can be `DEL` or `FEP`, defines the multipolygons geometry for the wildfire delineation for an AOI YY of the activation XXX.
- ``EMSRXXX_AOIYY_GRA_PRODUCT_naturalLandUseA.json`` defines the various multipolygons geometry for the grading damage levels for an AOI YY of the
activation XXX.

Those are the different grading levels of damage used in the Copernicus products:

Damage Level | CopernicusEMS class | EMS-98 class | Color
--- | --- |  --- | ---
0 | `No visible damage` | `No visible damage` | ![#000000](https://placehold.co/15x15/000000/000000.png)
1 | - | `Negligible to slight damage`  | ![#b4fe8e](https://placehold.co/15x15/b4fe8e/b4fe8e.png)
2 | `Possibly damaged` | `Moderate damage` | ![#fed98e](https://placehold.co/15x15/fed98e/fed98e.png)
3 | `Damaged` | `High damage` | ![#fe9929](https://placehold.co/15x15/fe9929/fe9929.png)
4 | `Destroyed` | `Destruction` | ![#cc4c02](https://placehold.co/15x15/cc4c02/cc4c02.png)



<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_DEL.png" width=30% height=30% alt>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_GRA.png" width=30% height=30% alt>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_S2L2A.png" width=30% height=30% alt>
    <br>
    <em>Here in order are reported the DEL map, GRA map and the actual sentinel-2 Image for the activation EMSR382 </em>
</p>


## Cloud masks

Creating cloud masks before making inferences on Sentinel-2 images is important because clouds can obscure or distort the underlying land cover or land use information that is the focus of the analysis. This can lead to inaccurate or incomplete results. Sentinel-2 images are often used for remote sensing applications, such as monitoring vegetation health, mapping land cover and land use, and detecting changes over time. However, clouds can interfere with these applications by blocking or reflecting the light that is captured by the satellite, which can result in missing or distorted data. By default, all images are retrieved from sentinel-hub with the condition of no more than 10 percent of cloud coverage. However some images have a relevant cloud coverage.

This dataset makes available for each image a cloud masks: the areas that are affected by clouds can be identified and excluded from future analyses. This ensures that the inferences made from the Sentinel-2 data are based on accurate and reliable information. The masks were generated using the CloudSen12 model.

The output prediction of CloudSen12 has 4 different layers for cloud coverage:

Label | Class | Class definitions | Color
--- | --- |  --- | ---
0 | `Clear` | Areas where no clouds are found | ![#67BEE0](https://placehold.co/15x15/67BEE0/67BEE0.png)
1 | `Clouds` | Clouds and heavy clouds are present, terrain is obscured and not visible  | ![#DCDCDC](https://placehold.co/15x15/DCDCDC/DCDCDC.png)
2 | `Light clouds` | Areas affected by light clouds, where could cause some issue in the terrain visibility. In this class are included also fog and wildfire's smoke.  | ![#B4B4B4](https://placehold.co/15x15/B4B4B4/B4B4B4.png)
3 | `Shadow` | This areas are in the shadow of the clouds. The terrain is partially/fully visible but the color and some bands of sentinel2 could be changed from real value. | ![#3C3C3C](https://placehold.co/15x15/3C3C3C/3C3C3C.png)

<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/images/EMSR382_AOI01_01_S2L2A.png" width=30% height=30% alt>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/images/EMSR382_AOI01_01_CM.png" width=30% height=30% alt>
    <br>
    <em>The cloud masks resulted from CloudSen12 model for activation EMSR382</em>
</p>

#### Smoothing-blend image patches

One challenge of using a U-Net for image segmentation is to have smooth predictions, especially if the receptive field of the neural network is a small amount of pixels. In the context of the U-Net architecture for image segmentation, blending image patches can be used to generate smooth predictions by reducing the effect of discontinuities at patch boundaries. This approach involves dividing the input image into overlapping patches, running the U-Net architecture on each patch individually, and then blending the resulting predictions together to form a single output image.

By blending the predictions from multiple patches, the resulting output image is typically smoother and more continuous than if a single U-Net model was trained on the entire input image. This can help to reduce artifacts and improve the overall quality of the segmentation results.
In this work the source code of cloudSen12 has been customized so that it could be smoothly predicted. The source code for smooth-blend is available [here](https://github.com/Vooban/Smoothly-Blend-Image-Patches)

<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/images/EMSR638_AOI01_01_S2L2A.png" width=30% height=30% alt>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/images/EMSR638_AOI01_01_CM.png" width=30% height=30% alt>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/images/EMSR638_AOI01_01_CMsmooth.png" width=30% height=30% alt>
    <br>
    <em>Cloud mask of activation EMSR638. It is clear that on left mask there are some problem due border effect in cloudSen12 model. On the right the result using smoothing</em>
</p>


## Landcover masks

In addition to wildfire delineation, severity and cloud masks, also the landcovers is provided for each image. In particular the models considered are: 
- [**ESRI 10m Annual Land Use Land Cover (2017-2021)**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1);
- [**ESRI 2020 Global Land Use Land Cover**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1);
- [**ESA WorldCover 10 m 2020**](https://esa-worldcover.org/en/data-access).

All this lancover are downloaded from [Planetary Computer](https://planetarycomputer.microsoft.com/)

All landcover models are based on sentinel2 10-meter resolution.

#### ESRI 10m Annual Land Use Land Cover (2017-2021) (9 classes)
There are 9 classes:

Original label | Remapped label | Class | Class definitions | Color
--- | --- | --- | --- | ---
0 | 255 | `No Data` | No data land | ![#FFFFFF](https://placehold.co/15x15/FFFFFF/FFFFFF.png)
1 | 0 | `Water` | Areas where water was predominantly present throughout the year; may not cover areas with sporadic or ephemeral water; contains little to no sparse vegetation, no rock outcrop nor built up features like docks; examples: rivers, ponds, lakes, oceans, flooded salt plains. | ![#1A5BAB](https://placehold.co/15x15/1A5BAB/1A5BAB.png)
2 | 1 | `Trees` | Any significant clustering of tall (~15 feet or higher) dense vegetation, typically with a closed or dense canopy; examples: wooded vegetation, clusters of dense tall vegetation within savannas, plantations, swamp or mangroves (dense/tall vegetation with ephemeral water or canopy too thick to detect water underneath). | ![#358221](https://placehold.co/15x15/358221/358221.png)
4 | 3 | `Flooded Vegetation` | Areas of any type of vegetation with obvious intermixing of water throughout a majority of the year; seasonally flooded area that is a mix of grass/shrub/trees/bare ground; examples: flooded mangroves, emergent vegetation, rice paddies and other heavily irrigated and inundated agriculture. | ![#87D19E](https://placehold.co/15x15/87D19E/87D19E.png)
5 | 4 | `Crops` | Human planted/plotted cereals, grasses, and crops not at tree height; examples: corn, wheat, soy, fallow plots of structured land. | ![#FFDB5C](https://placehold.co/15x15/FFDB5C/FFDB5C.png)
7 | 6 | `Built Area` | Human-made structures; major road and rail networks; large homogeneous impervious surfaces including parking structures, office buildings and residential housing; examples: houses, dense villages / towns / cities, paved roads, asphalt. | ![#ED022A](https://placehold.co/15x15/ED022A/ED022A.png)
8 | 7 | `Bare Ground` | Areas of rock or soil with very sparse to no vegetation for the entire year; large areas of sand and deserts with no to little vegetation; examples: exposed rock or soil, desert and sand dunes, dry salt flats/pans, dried lake beds, mines. | ![#EDE9E4](https://placehold.co/15x15/EDE9E4/EDE9E4.png)
9 | 8 | `Snow/Ice` | Large homogeneous areas of permanent snow or ice, typically only in mountain areas or highest latitudes; examples: glaciers, permanent snowpack, snow fields. | ![#F2FAFF](https://placehold.co/15x15/F2FAFF/F2FAFF.png)
10 | 9 | `Clouds` | No land cover information due to persistent cloud cover. | ![#C8C8C8](https://placehold.co/15x15/C8C8C8/C8C8C8.png)
11 | 10 | `Rangeland` | Open areas covered in homogeneous grasses with little to no taller vegetation; wild cereals and grasses with no obvious human plotting (i.e., not a plotted field); examples: natural meadows and fields with sparse to no tree cover, open savanna with few to no trees, parks/golf courses/lawns, pastures. | ![#C6AD8D](https://placehold.co/15x15/C6AD8D/C6AD8D.png)

<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_Annual9_LC.png" width=30% height=30% alt>
    <br>
    <em>Landcover ESRI 9 classes annual land use activation EMSR382_AOI01</em>
</p>

#### ESRI 2020 Global Land Use Land Cover
There are 10 classes:

Original label | Remapped label | Class | Class definitions | Color
--- | --- | --- | --- | ---
0 | 255 | `No Data` | No data land | ![#FFFFFF](https://placehold.co/15x15/FFFFFF/FFFFFF.png)
1 | 0 | `Water` | Water Areas where water was predominantly present throughout the year; may not cover areas with sporadic or ephemeral water; contains little to no sparse vegetation, no rock outcrop nor built up features like docks; examples: rivers, ponds, lakes, oceans, flooded salt plains. | ![#1A5BAB](https://placehold.co/15x15/1A5BAB/1A5BAB.png)
2 | 1 | `Trees` | Any significant clustering of tall (~15 feet or higher) dense vegetation, typically with a closed or dense canopy; examples: wooded vegetation, clusters of dense tall vegetation within savannas, plantations, swamp or mangroves (dense/tall vegetation with ephemeral water or canopy too thick to detect water underneath). | ![#358221](https://placehold.co/15x15/358221/358221.png)
3 | 2 | `Grass` | Open areas covered in homogeneous grasses with little to no taller vegetation; wild cereals and grasses with no obvious human plotting (i.e., not a plotted field); examples: natural meadows and fields with sparse to no tree cover, open savanna with few to no trees, parks/golf courses/lawns, pastures. | ![#A7D282](https://placehold.co/15x15/A7D282/A7D282.png)
4 | 3 | `Flooded Vegetation` | Areas of any type of vegetation with obvious intermixing of water throughout a majority of the year; seasonally flooded area that is a mix of grass/shrub/trees/bare ground; examples: flooded mangroves, emergent vegetation, rice paddies and other heavily irrigated and inundated agriculture. | ![#87D19E](https://placehold.co/15x15/87D19E/87D19E.png)
5 | 4 | `Crops` | Human planted/plotted cereals, grasses, and crops not at tree height; examples: corn, wheat, soy, fallow plots of structured land. | ![#FFDB5C](https://placehold.co/15x15/FFDB5C/FFDB5C.png)
6 | 5 | `Scrub/shrub` | Mix of small clusters of plants or single plants dispersed on a landscape that shows exposed soil or rock; scrub-filled clearings within dense forests that are clearly not taller than trees; examples: moderate to sparse cover of bushes, shrubs and tufts of grass, savannas with very sparse grasses, trees or other plants. | ![#EECFA8](https://placehold.co/15x15/EECFA8/EECFA8.png)
7 | 6 | `Built Area` |Human made structures; major road and rail networks; large homogeneous impervious surfaces including parking structures, office buildings and residential housing; examples: houses, dense villages / towns / cities, paved roads, asphalt. | ![#ED022A](https://placehold.co/15x15/ED022A/ED022A.png)
8 | 7 | `Bare Ground` | Areas of rock or soil with very sparse to no vegetation for the entire year; large areas of sand and deserts with no to little vegetation; examples: exposed rock or soil, desert and sand dunes, dry salt flats/pans, dried lake beds, mines. | ![#EDE9E4](https://placehold.co/15x15/EDE9E4/EDE9E4.png)
9 | 8 | `Snow/Ice` | Large homogeneous areas of permanent snow or ice, typically only in mountain areas or highest latitudes; examples: glaciers, permanent snowpack, snow fields. | ![#F2FAFF](https://placehold.co/15x15/F2FAFF/F2FAFF.png)
10 | 9 | `Clouds` | No land cover information due to persistent cloud cover. | ![#C8C8C8](https://placehold.co/15x15/C8C8C8/C8C8C8.png)

<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_Esri10_LC.png" width=30% height=30% alt>
    <br>
    <em>Landcover ESRI 10 classes 2020 land use activation EMSR382_AOI01</em>
</p>

#### ESA WorldCover 10m 2020
These are 10 classes:

Original label | Remapped label | Class | Class definitions | Color
--- | --- | --- | --- | ---
0 | 255 | `No Data` | No data land | ![#FFFFFF](https://placehold.co/15x15/FFFFFF/FFFFFF.png)
10 | 0 | `Trees` | This class includes any geographic area dominated by trees with a cover of 10% or more. Other land cover classes (shrubs and/or herbs in the understorey, built-up, permanent water bodies, …) can be present below the canopy, even with a density higher than trees. Areas planted with trees for afforestation purposes and plantations (e.g. oil palm, olive trees) are included in this class. This class also includes tree covered areas seasonally or permanently flooded with fresh water except for mangroves. | ![#006400](https://placehold.co/15x15/006400/006400.png)
20 | 1 | `Shrubland` | This class includes any geographic area dominated by natural shrubs having a cover of 10% or more. Shrubs are defined as woody perennial plants with persistent and woody stems and without any defined main stem being less than 5 m tall. Trees can be present in scattered form if their cover is less than 10%. Herbaceous plants can also be present at any density. The shrub foliage can be either evergreen or deciduous. | ![#FFBB22](https://placehold.co/15x15/FFBB22/FFBB22.png)
30 | 2 | `Grassland` | This class includes any geographic area dominated by natural herbaceous plants (Plants without persistent stem or shoots above ground and lacking definite firm structure): (grasslands, prairies, steppes, savannahs, pastures) with a cover of 10% or more, irrespective of different human and/or animal activities, such as: grazing, selective fire management etc. Woody plants (trees and/or shrubs) can be present assuming their cover is less than 10%. It may also contain uncultivated cropland areas (without harvest/ bare soil period) in the reference year | ![#FFFF4C](https://placehold.co/15x15/FFFF4C/FFFF4C.png)
40 | 3 | `Cropland` | Land covered with annual cropland that is sowed/planted and harvestable at least once within the 12 months after the sowing/planting date. The annual cropland produces an herbaceous cover and is sometimes combined with some tree or woody vegetation. Note that perennial woody crops will be classified as the appropriate tree cover or shrub land cover type. Greenhouses are considered as built-up. | ![#F096FF](https://placehold.co/15x15/F096FF/F096FF.png)
50 | 4 | `Built-up` | Land covered by buildings, roads and other man-made structures such as railroads. Buildings include both residential and industrial building. Urban green (parks, sport facilities) is not included in this class. Waste dump deposits and extraction sites are considered as bare. | ![#FF0000](https://placehold.co/15x15/FF0000/FF0000.png)
60 | 5 | `Bare/sparse vegetation` | Lands with exposed soil, sand, or rocks and never has more than 10% vegetated cover during any time of the year | ![#B4B4B4](https://placehold.co/15x15/B4B4B4/B4B4B4.png)
70 | 6 | `Snow and Ice` |This class includes any geographic area covered by snow or glaciers persistently | ![#F0F0F0](https://placehold.co/15x15/F0F0F0/F0F0F0.png)
80 | 7 | `Permanent water bodies` | This class includes any geographic area covered for most of the year (more than 9 months) by water bodies: lakes, reservoirs, and rivers. Can be either fresh or salt-water bodies. In some cases the water can be frozen for part of the year (less than 9 months). | ![#0064C8](https://placehold.co/15x15/0064C8/0064C8.png)
90 | 8 | `Herbaceous wetland` | Land dominated by natural herbaceous vegetation (cover of 10% or more) that is permanently or regularly flooded by fresh, brackish or salt water. It excludes unvegetated sediment (see 60), swamp forests (classified as tree cover) and mangroves see 95) | ![#0096A0](https://placehold.co/15x15/0096A0/0096A0.png)
95 | 9 | `Mangroves` | Taxonomically diverse, salt-tolerant tree and other plant species which thrive in intertidal zones of sheltered tropical shores, "overwash" islands, and estuaries. | ![#00CF75](https://placehold.co/15x15/00CF75/00CF75.png)
100 | 10 | `Moss and lichen` | Land covered with lichens and/or mosses. Lichens are composite organisms formed from the symbiotic association of fungi and algae. Mosses contain photo-autotrophic land plants without true leaves, stems, roots but with leaf-and stemlike organs. | ![#FAE6A0](https://placehold.co/15x15/FAE6A0/FAE6A0.png)

More informations are availble in the [Esa worldcover manual](https://esa-worldcover.s3.amazonaws.com/v100/2020/docs/WorldCover_PUM_V1.0.pdf)

<p align="center">
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_ESA_LC.png" width=30% height=30% alt>
    <br>
    <em>Landcover ESA worldcover 2020 land use activation EMSR382_AOI01</em>
</p>

