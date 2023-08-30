# Python file with general settings

# DEBUG
DEBUGDownload = False

# Sentinel client info
# Register to https://apps.sentinel-hub.com/dashboard/#/ --> User Setting --> OAuth clients --> Create new client name
sh_client_id = "df6dc758-639c-4bee-9457-25b454d53974"
sh_client_secret = "p4I?!2GcF7Mw_+)/,&d#^Rp^Pou,<;N;:tWDaJsq"
planetary_computer_API_KEY = "1ded7aabaa5b4fe6b69f26fadbe0ae7b"

# Retrivial 
resolution = 10             # 10 meter resolution
cloudCoverMax = 0.1        # max cloud coverage tollerance 10 percent
min_dims = 512              # min dimension tile download from sentinelHub
splitMaxDimension = 2028    # max dimension tile download from sentinelHub

# Grading                            Copernicus EMS classes     EMS-98 classes
damageLevel = {"destruction" :          ["Destroyed" ,          "Destruction"],                        #level 4
                "high damage" :         ["Damaged",             "High damage"],                        #level 3
                "moderate damage" :     ["Possibly damaged",    "Moderate damage"],                    #level 2
                "negligible damage":    [                       "Negligible to slight damage"],        #level 1
                "not analized":         ["Not analyzed",        "No visible damage"]                   #level 0
            } 

# Path
copernicusFolderPath =  "copernicusData/cems-rm-viewer/" #"copernicusData/cems-rm-viewer/" # "copernicusData/colomba/" #"copernicusData/cems-rm-viewer/"
root_folder = "../../"
output_folder = root_folder + "data/" #"output/" #"data/"
RGB_folder = root_folder + "dataRGB/"
imagesFolder = "assets/images/"

dateJson = "date.json"
evalscript = "evalScript.js"

datasetPreConfigured = "csv_files/datasetPreConfigured.csv" # "csv_files/datasetPreConfiguredColomba.csv"
satelliteData = "csv_files/satelliteData.csv"
cloudCoveragecsv = "csv_files/cloudCoverage.csv"
logPath = "csv_files/log.txt"       

# Custom cloud percentage
customCloudPercentage = False
cloudPercentage = 6500          # level 1 heavy clouds
lightCloudPercentage = 7000     # level 2 light clouds

# Image RGB
contrast_coeff = 10000/255*3.5*5.7     # image color contrast 
brightness_coeff = 30 #20              # image brightness 

# Flag images
saveMergedTif = False       # save a single merged tiff per AoI ( it can take several GB of memory )
saveRGBImage = True         # save tif image also in .png
saveOverlapMask = False     # save overlap mask between burned area and cloud
saveGradingRGB = False      # save different level of grading in separate mask as png
saveCloudImage = False      # save cloud tif prediction from cluodSen12

# LandCover
landCover =  ["ESA", "LC9", "LC10"] #"ESA" #"LC10" #"LC9"

