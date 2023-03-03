# Python file with general variables

# DEBUG
DEBUGDownload = False

# Sentinel client info
sh_client_id = "96b24011-6266-400e-9121-feb214296ebb"
sh_client_secret = "b2{uHX,me%kM6PK3ix<*9W[KPbAQbt.c.rB#sP&;"

# Retrivial 
min_dims = 512
resolution = 10             # 10 meter resolution
cloudCoverMax = 0.1         # max cloud coverage tollerance 10 percent
splitMaxDimension = 2028    # max dimension tile download

# Grading
damageLevel = {"destruction" :          ["Destroyed" , "Destruction"],                              #level 4
                "high damage" :         ["Damaged", "High damage"],                                 #level 3
                "moderate damage" :     ["Possibly damaged", "Moderate damage"],                    #level 2
                "negligible damage":    ["Negligible to slight damage"],                            #level 1
                "not analized":         ["Not analyzed", "No visible damage"]                       #level 0
            }

# Path
copernicusFolderPath = "copernicusData/cems-rm-viewer/"
root_folder = "../../"
output_folder = root_folder + "data/" #"output/" #"data/"
RGB_folder = root_folder + "dataRGB/"
imagesFolder = "assets/images/"

dateJson = "date.json"
evalscript = "evalScript.js"

datasetPreConfigured = "csv_files/datasetPreConfigured.csv"
satelliteData = "csv_files/satelliteData.csv"
cloudCoveragecsv = "csv_files/cloudCoverage.csv"
logPath = "csv_files/log.txt"       

# Cloud percentage
customCloudPercentage = False
cloudPercentage = 6500          # level 1 heavy clouds
lightCloudPercentage = 7000     # level 2 light clouds

# Image RGB
contrast_coeff = 10000/255*3.5*5.5      # image color contrast 
brightness_coeff = -10 #-20             # image brightness 

# Flag images
saveMergedTif = False       # save a single merged tiff per AoI ( it can take several GB of memory )
saveRGBImage = True         # save tif image also in .png
saveOverlapMask = False     # save overlap mask between burned area and cloud
saveGradingRGB = False      # save different level of grading in separate mask as png
saveCloudImage = False      # save cloud tif prediction from cluodSen12

# LandCover
landCover =  ["ESA", "LC9", "LC10"] #"ESA" #"LC10" #"LC9"

