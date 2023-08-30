import cloudsen12
import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
import xarray as xr
import gc

import torch
import scipy.signal
from tqdm import tqdm

from maskay.tensorsat import TensorSat
from maskay.predict import Predictor
from maskay.module import MaskayModule

from maskay.utils import softmax
from maskay.library.unetmobv2 import model_setup


class CloudCover:

    def __init__(self, cropsize: int, overlap: int):
        self.cloudBands = ["Aerosol", "Blue", "Green", "Red", "RedEdge1", "RedEdge2", "RedEdge3", "NIR", "NIR2", "WaterVapor", "SWIR1", "SWIR2"]
        
        # CloudSen12 with Smoothing
        self.cloud_model = UnetMobV2Custom()
        self.cloud_predictor = PredictorCustom(
            cropsize = cropsize,
            overlap = overlap,
            device = "cpu", #"cpu" #"cuda"
            quiet = True
        )
        
        # Uncomment for CloudSen12 without Smoothing
        # self.cloud_model = cloudsen12.UnetMobV2()
        # self.cloud_predictor = cloudsen12.Predictor(
        #     cropsize = cropsize,
        #     overlap = overlap,
        #     device = "cpu",
        #     quiet = True
        # )


    def predictCover(self, image: np.ndarray):

        """Predicts the cloud cover in a given satellite image.
        
        Args: 
            image_mask (np.ndarray):A numpy ndarray representing the input satellite image.
            Dimensions are expected to be [number_of_bands, height, width].
        
        Returns: 
            np.ndarray: a prediction of cloud cover for the provided image."""

        cirrus = np.full((image.shape[1], image.shape[2]), 0.0001) # not zero value
        image = np.nan_to_num(image)
        
        cloud_bands_dict = {"Cirrus" : cirrus} # create 13th dummy layer 
        for i, band in enumerate(self.cloudBands):
            cloud_bands_dict[band] = image[i] 
        
        # create TensorSat dict: key -> bandName, value -> xr.DataArray
        cloud_tensor = cloudsen12.TensorSat(**cloud_bands_dict, cache=True, align=False)

        prediction = self.cloud_predictor.predictCustom(self.cloud_model, cloud_tensor)     # with smoothing
        # prediction = self.cloud_predictor.predict(self.cloud_model, cloud_tensor)         # without smoothing
        return prediction


#########################################################################################################
                                                # MASKAY PACKAGE
#########################################################################################################
# this is the customization for namespace package Maskay
# original code can be found here: https://github.com/cloudsen12/maskay/tree/main/maskay
class MaskayModuleCustom(MaskayModule):

    def __init__(self, cropsize, overlap, device, batchsize, order, quiet):
        super().__init__(cropsize, overlap, device, batchsize, order, quiet)
        self.cropsize = cropsize
        self.overlap = overlap
        self.device = device
        self.batchsize = batchsize
        self.order = order
        self.quiet = quiet


    def _predictCustom(self, tensor: TensorSat):
        
        input_nb_classes = 13
        output_nb_classes = 4 # classes: 0 = clean sky, 1 = clouds, 2 = light_clouds/fog/smoke, 3 = cloud's shadow
        window_size = 512
        smooth = SmoothImage()

        # Raster ref (lowest resolution)
        rbase = tensor.rasterbase()

        # ['Aerosol', 'Blue', 'Green', 'Red', 'RedEdge1', 'RedEdge2', 'RedEdge3', 'NIR', 'NIR2', 'WaterVapor', 'Cirrus', 'SWIR1', 'SWIR2']
        level_names = tensor.dictarray.keys()

        tensorList = list() 

        for key in level_names:
            tensorXarrayImg = tensor.dictarray[key]
            tensorList.append(tensorXarrayImg.to_numpy())

        # create the 13 bands image
        tensorImage = np.stack(tensorList, axis=0) # (13, H, W)
        
        # channel last
        tensorImage = np.transpose(tensorImage, (1, 2, 0)) # (H, W, 13)

        # preprocessing of the image
        tensorImage = self.inProcessing(tensorImage)

        # smoothing of the image with spline interpolation
        outTensorImage = smooth.predict_img_with_smooth_windowing(tensorImage,
                                                window_size = window_size,
                                                subdivisions = 2,  # Minimal amount of overlap for windowing. Must be an even number.
                                                input_nb_classes = input_nb_classes,
                                                output_nb_classes = output_nb_classes,
                                                pred_func = (
                                                    lambda img_batch_subdiv: self._run(img_batch_subdiv))
                                                )
        
        # shape outTensorImage: (H, W, 4)
        # channel first
        outTensorImage = np.transpose(outTensorImage, (2, 0, 1)) # (4, H, W)

        outTensorImage = np.expand_dims(outTensorImage, 0) # (1, 4, H, W)

        # postprocessing of the image
        outTensorImage = self.outProcessing(outTensorImage)

        outTensorImage = np.squeeze(outTensorImage) # (4, H, W))

        # Create the output tensor rio xarray object
        xcoord, ycoord = rbase.coords["y"].values, rbase.coords["x"].values
        coords = [np.arange(0, output_nb_classes), xcoord, ycoord]
        dims = ["band", "y", "x"]

        return (
            xr.DataArray(outTensorImage, coords=coords, dims=dims)
            .rio.write_nodata(-999)
            .rio.write_transform(rbase.rio.transform())
        )
 

class ModuleCustom(MaskayModuleCustom):

    def __init__(
        self,
        cropsize: int = 512,
        overlap: int = 32,
        device: str = "cpu",
        batchsize: int = 1,
        order: str = "BCHW",
        quiet: int = False,
    ):
        super().__init__(cropsize, overlap, device, batchsize, order, quiet)

    def inProcessing(self, tensor: np.ndarray):
        pass

    def forward(self, x):
        pass

    def outProcessing(self, tensor: np.ndarray):
        pass

    def _run(self, tensor: np.ndarray):
        tensor = torch.Tensor(tensor)

        # Run the model
        with torch.no_grad():
            if self.device == "cuda":
                tensor = tensor.cuda()
            tensor = self.forward(tensor).detach().cpu().numpy()
            torch.cuda.empty_cache()
        return tensor


class UnetMobV2Custom(ModuleCustom):
    def __init__(self):
        super().__init__()
        self.model = model_setup()

    def forward(self, x):
        return self.model(x)

    def inProcessing(self, tensor: np.ndarray):
        # If all the pixels are zero skip the run and outProcessing.
        if np.sum(tensor) == 0:
            shp = tensor.shape
            tensor = np.zeros(
                (shp[0], 4, shp[2], shp[3])
            ) # 4 is the number of the output classes
            return [tensor]
        return tensor / 10000

    def outProcessing(self, tensor: np.ndarray):
        return (softmax(tensor, axis=1) * 10000).astype(np.int16)


class PredictorCustom(Predictor):
    def __init__(
        self,
        cropsize: int = 512,
        overlap: int = 32,
        device: str = "cpu",
        batchsize: int = 1,
        order: int = "BCHW",
        quiet: int = False,
    ):
        super().__init__(cropsize, overlap, device, batchsize, order, quiet)
        self.batchsize = batchsize
        self.cropsize = cropsize
        self.overlap = overlap
        self.device = device
        self.order = order
        self.quiet = quiet
        self.result = None

    def predictCustom(self, model, tensor: TensorSat):
        model.batchsize = self.batchsize
        model.cropsize = self.cropsize
        model.overlap = self.overlap
        model.device = self.device
        model.order = self.order
        model.quiet = self.quiet
        self.result = model._predictCustom(tensor=tensor)
        return self.result


#########################################################################################################
                                                # SMOOTHING 
#########################################################################################################
# Source code here: https://github.com/Vooban/Smoothly-Blend-Image-Patches/blob/master/smooth_tiled_predictions.py
# Video explanation: https://www.youtube.com/watch?v=HrGn4uFrMOM&t=3s&ab_channel=DigitalSreeni
# Wiki: https://medium.com/vooban-ai/satellite-image-segmentation-a-workflow-with-u-net-7ff992b2a56e
PLOT_PROGRESS = True
cached_2d_windows = dict()
class SmoothImage:
    def __init__(self) -> None:
        pass

    def _spline_window(self, window_size: int, power: int = 2):

        """Generates a 1-dimensional spline of order 'power' (typically 2), in the designated window.
        https://www.wolframalpha.com/input/?i=y%3Dx**2,+y%3D-(x-2)**2+%2B2,+y%3D(x-4)**2,+from+y+%3D+0+to+2
        
        Args:
            window_size (int): size of the interested window power (int, optional): Order of the spline. Defaults to 2.
        
        Returns:
            np.ndarray: 1D spline
        """

        intersection = int(window_size/4)
        wind_outer = (abs(2*(scipy.signal.triang(window_size))) ** power)/2
        wind_outer[intersection:-intersection] = 0

        wind_inner = 1 - (abs(2*(scipy.signal.triang(window_size) - 1)) ** power)/2
        wind_inner[:intersection] = 0
        wind_inner[-intersection:] = 0

        wind = wind_inner + wind_outer
        wind = wind / np.average(wind)

        return wind

    def _window_2D(self, window_size: int, power: int = 2):

        """Makes a 1D window spline function, then combines it to return a 2D window function.
        The 2D window is useful to smoothly interpolate between patches.

        Args:
            window_size (int): size of the window (patch)
            power (int, optional): Which order for the spline. Defaults to 2.

        Returns:
            np.ndarray: numpy array containing a 2D spline function
        """

        # Memorization to avoid remaking it for every call
        # since the same window is needed multiple times
        global cached_2d_windows
        key = "{}_{}".format(window_size, power)
        if key in cached_2d_windows:
            wind = cached_2d_windows[key]
        else:
            wind = self._spline_window(window_size, power)
            wind = np.expand_dims(wind, 1)
            wind = np.expand_dims(wind, 1)
            wind = wind * wind.transpose(1, 0, 2)
 
            cached_2d_windows[key] = wind
        return wind


    def _pad_img(self, img: np.ndarray, window_size: int, subdivisions: int):

        """Add borders to the given image for a "valid" border pattern according to "window_size" and "subdivisions".
        Image is expected as a numpy array with shape (width, height, channels).

        Args:
            image (np.ndarray): input image, 3D channels-last tensor
            window_size (int): size of a single patch, useful to compute padding
            subdivisions (int): amount of overlap, useful for padding

        Returns:
            np.ndarray: same image, padded specularly by a certain amount in every direction
        """

        # compute the pad as (window - window/subdivisions)
        aug = int(round(window_size * (1 - 1.0/subdivisions)))
        more_borders = ((aug, aug), (aug, aug), (0, 0))
        ret = np.pad(img, pad_width=more_borders, mode='reflect')

        return ret


    def _unpad_img(self, padded_img: np.ndarray, window_size: int, subdivisions: int):

        """Reverts changes made by 'pad_image'. The same padding is removed, so window_size and subdivisions
        must be coherent.

        Args:
            padded_image (np.ndarray): image with padding still applied
            window_size (int): size of a single patch
            subdivisions (int): subdivisions to compute overlap

        Returns:
            np.ndarray: image without padding, 2D channels-last tensor
        """
        # compute the pad as (window - window/subdivisions)
        aug = int(round(window_size * (1 - 1.0/subdivisions)))
        ret = padded_img[
            aug:-aug,
            aug:-aug,
            :
        ]
        return ret


    def _rotate_mirror_do(self, im: np.ndarray):

        """Duplicates an image with shape (h, w, channels) 8 times, in order
        to have all the possible rotations and mirrors of that image that fits the
        possible 90 degrees rotations. https://en.wikipedia.org/wiki/Dihedral_group

        Args:
            image (np.ndarray): input image, already padded.

        Returns:
            List[np.ndarray]: list of images, rotated and mirrored.
        """

        mirrs = []
        mirrs.append(np.array(im))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=1))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=2))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=3))
        im = np.array(im)[:, ::-1]
        mirrs.append(np.array(im))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=1))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=2))
        mirrs.append(np.rot90(np.array(im), axes=(0, 1), k=3))
        return mirrs


    def _rotate_mirror_undo(self, im_mirrs: np.ndarray):

        """Reverts the 8 duplications provided by rotate and mirror.
        This restores the transformed inputs to the original position, then averages them.

        Args:
            variants (List[np.ndarray]): D4 dihedral group of the same image

        Returns:
            np.ndarray: averaged result over the given input.
        """
        origs = []
        origs.append(np.array(im_mirrs[0]))
        origs.append(np.rot90(np.array(im_mirrs[1]), axes=(0, 1), k=3))
        origs.append(np.rot90(np.array(im_mirrs[2]), axes=(0, 1), k=2))
        origs.append(np.rot90(np.array(im_mirrs[3]), axes=(0, 1), k=1))
        origs.append(np.array(im_mirrs[4])[:, ::-1])
        origs.append(np.rot90(np.array(im_mirrs[5]), axes=(0, 1), k=3)[:, ::-1])
        origs.append(np.rot90(np.array(im_mirrs[6]), axes=(0, 1), k=2)[:, ::-1])
        origs.append(np.rot90(np.array(im_mirrs[7]), axes=(0, 1), k=1)[:, ::-1])
        return np.mean(origs, axis=0)


    def _windowed_subdivs(self, padded_img: np.ndarray, window_size: int, subdivisions: int, output_nb_classes: int, pred_func):
        
        """Create tiled overlapping patches.

        
        Note:
            patches_resolution_along_X == patches_resolution_along_Y == window_size

        Args:
            padded_image (np.ndarray): image with padding still applied
            window_size (int): size of a single patch
            subdivisions (int): amount of overlap, useful for padding
            output_nb_classes (int): number output classes 
            pred_func: function object reference

        Returns:
            np.ndarray: 5D numpy array of shape = (
                nb_patches_along_X,
                nb_patches_along_Y,
                patches_resolution_along_X,
                patches_resolution_along_Y,
                nb_output_channels
            )
        """

        WINDOW_SPLINE_2D = self._window_2D(window_size=window_size, power=2)

        step = int(window_size/subdivisions)
        padx_len = padded_img.shape[0]
        pady_len = padded_img.shape[1]
        subdivsList = []

        count = 0
        for i in range(0, padx_len-window_size+1, step):
            subdivsList.append([])
            for j in range(0, pady_len-window_size+1, step):
                patch = padded_img[i:i+window_size, j:j+window_size, :]
                subdivsList[-1].append(patch)
                count += 1

        
        # Here, `gc.collect()` clears RAM between operations.
        # It should run faster if they are removed, if enough memory is available.
        # NOTE: the first two dimensions of the subdivision represent the height and width of the padded image considering a tile of 512x512 (bear in mind the overlap)
        gc.collect()
        subdivs = np.array(subdivsList) # (5, 4, 512, 512, 13)

        gc.collect()
        a, b, c, d, e = subdivs.shape

        subdivs = subdivs.reshape(a * b, c, d, e) # (20, 512, 512, 13)
        gc.collect()

        subdivs = np.transpose(subdivs, (0, 3, 1, 2)) # (20, 13, 512, 512)

        # the cloudSen12 model want a channel first 
        subdivs = pred_func(subdivs) # Here cloudSen12 prediction
        
        subdivs = np.transpose(subdivs, (0, 2, 3, 1)) # (20, 512, 512, 13)

        gc.collect()
        subdivs = np.array([patch * WINDOW_SPLINE_2D for patch in subdivs]) # (20, 512, 512, 4)

        gc.collect()

        # Such 5D array:
        subdivs = subdivs.reshape(a, b, c, d, output_nb_classes) # (5, 4, 512, 512, 4)
        gc.collect()

        return subdivs


    def _recreate_from_subdivs(self, subdivs: np.ndarray, window_size: int, subdivisions: int, padded_out_shape):
        
        """Merge tiled overlapping patches smoothly.

        Args:
            subdivs (np.ndarray): image with padding still applied
            window_size (int): size of a single patch
            subdivisions (int): amount of overlap, useful for padding
            padded_out_shape (List): shape of the padded image ( channel first )

        Returns:
            np.ndarray: one padded result
        """

        step = int(window_size/subdivisions)
        padx_len = padded_out_shape[0]
        pady_len = padded_out_shape[1]

        y = np.zeros(padded_out_shape)

        a = 0
        for i in range(0, padx_len-window_size+1, step):
            b = 0
            for j in range(0, pady_len-window_size+1, step):
                windowed_patch = subdivs[a, b]
                y[i:i+window_size, j:j+window_size] = y[i:i+window_size, j:j+window_size] + windowed_patch
                b += 1
            a += 1
        return y / (subdivisions ** 2)


    def predict_img_with_smooth_windowing(self, input_img, window_size, subdivisions, input_nb_classes, output_nb_classes, pred_func):
        
        """Apply the `pred_func` function to square patches of the image, and overlap
        the predictions to merge them smoothly.
        See 6th, 7th and 8th idea here:
        http://blog.kaggle.com/2017/05/09/dstl-satellite-imagery-competition-3rd-place-winners-interview-vladimir-sergey/
        """

        pad = self._pad_img(input_img, window_size, subdivisions)
        pads = self._rotate_mirror_do(pad)

        # Note that the implementation could be more memory-efficient by merging
        # the behavior of `_windowed_subdivs` and `_recreate_from_subdivs` into
        # one loop doing in-place assignments to the new image matrix, rather than
        # using a temporary 5D array.

        # It would also be possible to allow different (and impure) window functions
        # that might not tile well. Adding their weighting to another matrix could
        # be done to later normalize the predictions correctly by dividing the whole
        # reconstructed thing by this matrix of weightings - to normalize things
        # back from an impure windowing function that would have badly weighted
        # windows.

        # For example, since the U-net of Kaggle's DSTL satellite imagery feature
        # prediction challenge's 3rd place winners use a different window size for
        # the input and output of the neural net's patches predictions, it would be
        # possible to fake a full-size window which would in fact just have a narrow
        # non-zero dommain. This may require to augment the `subdivisions` argument
        # to 4 rather than 2.

        res = []
        for pad in tqdm(pads):
            # For every rotation:
            sd = self._windowed_subdivs(pad, window_size, subdivisions, output_nb_classes, pred_func)
            one_padded_result = self._recreate_from_subdivs(
                sd, window_size, subdivisions,
                padded_out_shape=list(pad.shape[:-1])+[output_nb_classes]) 

            res.append(one_padded_result)

        # Merge after rotations:
        padded_results = self._rotate_mirror_undo(res)

        prediction = self._unpad_img(padded_results, window_size, subdivisions)

        prediction = prediction[:input_img.shape[0], :input_img.shape[1], :]

        return prediction


 