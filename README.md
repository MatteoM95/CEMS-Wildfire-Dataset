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

In addition to wildfire delineation mask and clod masks, also landcovers for each image are create. In particular the model considered are: 
- [**ESRI 10m Annual Land Use Land Cover (2017-2021)**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1);
- [**ESRI 2020 Global Land Use Land Cover**](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9553499&tag=1);
- [**ESA WorldCover 10 m 2020**](https://esa-worldcover.org/en/data-access).
All landcover models are based on sentinel2 10-meter resolution images.

#### ESRI 10m Annual Land Use Land Cover (2017-2021) (9 classes)
There are 9 classes:

Original label | Remapped label | Class | Class definitions | Color
--- | --- | --- | --- | ---
0 | 255 | `No Data` | No data land | ![#FFFFFF](https://placehold.co/15x15/FFFFFF/FFFFFF.png)
1 | 0 | `Water` | Water Areas where water was predominantly present throughout the year; may not cover areas with sporadic or ephemeral water; contains little to no sparse vegetation, no rock outcrop nor built up features like docks; examples: rivers, ponds, lakes, oceans, flooded salt plains. | ![#1A5BAB](https://placehold.co/15x15/1A5BAB/1A5BAB.png)
2 | 1 | `Trees` | Any significant clustering of tall (~15 feet or higher) dense vegetation, typically with a closed or dense canopy; examples: wooded vegetation, clusters of dense tall vegetation within savannas, plantations, swamp or mangroves (dense/tall vegetation with ephemeral water or canopy too thick to detect water underneath). | ![#358221](https://placehold.co/15x15/358221/358221.png)
4 | 3 | `Flooded Vegetation` | Areas of any type of vegetation with obvious intermixing of water throughout a majority of the year; seasonally flooded area that is a mix of grass/shrub/trees/bare ground; examples: flooded mangroves, emergent vegetation, rice paddies and other heavily irrigated and inundated agriculture. | ![#87D19E](https://placehold.co/15x15/87D19E/87D19E.png)
5 | 4 | `Crops` | Human planted/plotted cereals, grasses, and crops not at tree height; examples: corn, wheat, soy, fallow plots of structured land. | ![#FFDB5C](https://placehold.co/15x15/FFDB5C/FFDB5C.png)
7 | 6 | `Built Area` |Human made structures; major road and rail networks; large homogeneous impervious surfaces including parking structures, office buildings and residential housing; examples: houses, dense villages / towns / cities, paved roads, asphalt. | ![#ED022A](https://placehold.co/15x15/ED022A/ED022A.png)
8 | 7 | `Bare Ground` | Areas of rock or soil with very sparse to no vegetation for the entire year; large areas of sand and deserts with no to little vegetation; examples: exposed rock or soil, desert and sand dunes, dry salt flats/pans, dried lake beds, mines. | ![#EDE9E4](https://placehold.co/15x15/EDE9E4/EDE9E4.png)
9 | 8 | `Snow/Ice` | Large homogeneous areas of permanent snow or ice, typically only in mountain areas or highest latitudes; examples: glaciers, permanent snowpack, snow fields. | ![#F2FAFF](https://placehold.co/15x15/F2FAFF/F2FAFF.png)
10 | 9 | `Clouds` | No land cover information due to persistent cloud cover. | ![#C8C8C8](https://placehold.co/15x15/C8C8C8/C8C8C8.png)
11 | 10 | `Rangeland` | Open areas covered in homogeneous grasses with little to no taller vegetation; wild cereals and grasses with no obvious human plotting (i.e., not a plotted field); examples: natural meadows and fields with sparse to no tree cover, open savanna with few to no trees, parks/golf courses/lawns, pastures. | ![#C6AD8D](https://placehold.co/15x15/C6AD8D/C6AD8D.png)

<p>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_Annual9_LC.png" width=50% height=50% alt>
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

<p>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_Esri10_LC.png" width=50% height=50% alt>
    <em>Landcover ESRI 10 classes 2020 land use activation EMSR382_AOI01</em>
</p>

#### ESA WorldCover 10m 2020
These are 10 classes:

Original label | Remapped label | Class | Class definitions | Color
--- | --- | --- | --- | ---
0 | 255 | `No Data` | No data land | ![#FFFFFF](https://placehold.co/15x15/FFFFFF/FFFFFF.png)
10 | 0 | `Trees` | This class includes any geographic area dominated by trees with a cover of 10% or more. Other land cover classes (shrubs and/or herbs in the understorey, built-up, permanent water bodies, …) can be present below the canopy, even with a density higher than trees.
Areas planted with trees for afforestation purposes and plantations (e.g. oil palm, olive trees) are included in this class. This class also includes tree covered areas seasonally or permanently flooded with fresh water except for mangroves. | ![#006400](https://placehold.co/15x15/006400/006400.png)
20 | 1 | `Shrubland` | This class includes any geographic area dominated by natural shrubs having a cover of 10% or more. Shrubs are defined as woody perennial plants with persistent and woody stems and without any defined main stem being less than 5 m tall. Trees can be present in scattered form if their cover is less than 10%. Herbaceous plants can also be present at any density. The shrub foliage can be either evergreen or deciduous. | ![#FFBB22](https://placehold.co/15x15/FFBB22/FFBB22.png)
30 | 2 | `Grassland` | This class includes any geographic area dominated by natural herbaceous plants (Plants without persistent stem or shoots above ground and lacking definite firm structure): (grasslands, prairies, steppes, savannahs, pastures) with a cover of 10% or more, irrespective of different human and/or animal activities, such as: grazing, selective fire management etc. Woody plants (trees and/or shrubs) can be present assuming their cover is less than 10%. It may also contain uncultivated cropland areas (without harvest/ bare soil period) in the reference year | ![#FFFF4C](https://placehold.co/15x15/FFFF4C/FFFF4C.png)
40 | 3 | `Cropland` | Land covered with annual cropland that is sowed/planted and harvestable at least once within the 12 months after the sowing/planting date. The annual cropland produces an herbaceous cover and is sometimes combined with some tree or woody vegetation. Note that perennial woody crops will be classified as the appropriate tree cover or shrub land cover type. Greenhouses are considered as built-up. | ![#F096FF](https://placehold.co/15x15/F096FF/F096FF.png)
50 | 4 | `Built-up` | Land covered by buildings, roads and other man-made structures such as railroads. Buildings include both residential and industrial building. Urban green (parks, sport facilities) is not included in this class. Waste dump deposits and extraction sites are considered as bare. | ![#FF0000](https://placehold.co/15x15/FF0000/FF0000.png)
60 | 5 | `Bare/sparse vegetation` | Lands with exposed soil, sand, or rocks and never has more than 10% vegetated cover during any time of the year | ![#B4B4B4](https://placehold.co/15x15/B4B4B4/B4B4B4.png)
70 | 6 | `Snow and Ice` |This class includes any geographic area covered by snow or glaciers
persistently | ![#F0F0F0](https://placehold.co/15x15/F0F0F0/F0F0F0.png)
80 | 7 | `Permanent water bodies` | This class includes any geographic area covered for most of the year (more than 9 months) by water bodies: lakes, reservoirs, and rivers. Can be either fresh or salt-water bodies. In some cases the water can be frozen for part of the year (less than 9 months). | ![#0064C8](https://placehold.co/15x15/0064C8/0064C8.png)
90 | 8 | `Herbaceous wetland` | Land dominated by natural herbaceous vegetation (cover of 10% or
more) that is permanently or regularly flooded by fresh, brackish or salt water. It excludes unvegetated sediment (see 60), swamp forests (classified as tree cover) and mangroves see 95) | ![#0096A0](https://placehold.co/15x15/0096A0/0096A0.png)
95 | 9 | `Mangroves` | Taxonomically diverse, salt-tolerant tree and other plant species
which thrive in intertidal zones of sheltered tropical shores, "overwash" islands, and estuaries. | ![#00CF75](https://placehold.co/15x15/00CF75/00CF75.png)
100 | 10 | `Moss and lichen` | Land covered with lichens and/or mosses. Lichens are composite
organisms formed from the symbiotic association of fungi and algae. Mosses contain photo-autotrophic land plants without true leaves, stems, roots but with leaf-and stemlike organs. | ![#FAE6A0](https://placehold.co/15x15/FAE6A0/FAE6A0.png)

All informations are availble in the [Esa worldcover manual](https://esa-worldcover.s3.amazonaws.com/v100/2020/docs/WorldCover_PUM_V1.0.pdf)

<p>
    <img src="https://github.com/MatteoM95/CEMS-Wildfire-Dataset/blob/main/assets/sample/EMSR382/AOI01/EMSR382_AOI01_01/EMSR382_AOI01_01_ESA_LC.png" width=50% height=50% alt>
    <em>Landcover ESA worldcover 2020 land use activation EMSR382_AOI01</em>
</p>

in case broken git: git reset --hard origin/master
