# CEMS-Wildfire-Dataset
Copernicus emergency management service wildfire dataset from 2019 to 2022.

# Introduction
This dataset wiil be part of European project [SAFERS](https://safers-project.eu/) ‘Structured Approaches for Forest Fire Emergencies in Resilient Societies’. It is going to create an open and integrated platform featuring a forest fire Decision Support System. The platform will use information from different sources: earth observations from Copernicus and GEOSS, fire sensors in forests, topographic data, weather forecasts and even crowdsourced data from social media and other apps that can be used by citizens and first responders to provide situational in-field information.
A sample from dataset is available [here](ZZZZ)


# Installation packages
First install all requirements in [assets/requirements.txt](ZZZZZ) and cloudSen12 library
```
pip install -r /assets/requirements.txt
pip install cloudsen12
```
In case error with cartopy library, install
```
sudo apt -y install libgeos-dev
```


# How to run and create the dataset
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


# Landcovers

[Esa worldcover manual](https://esa-worldcover.s3.amazonaws.com/v100/2020/docs/WorldCover_PUM_V1.0.pdf)

in case broken git: git reset --hard origin/master