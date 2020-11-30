import numpy as np
import operator
from scipy.special import gamma
from scipy import ndimage
import pylab
import re
from os import path
from math import asin, tan, isnan
from image_tools import draw_sphere_mm, get_translation
from gate_macro_parser.phantom_parser import JPET
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from optparse import OptionParser


parser = OptionParser()

parser.add_option("-i", "--inter_file",
                  action="store", dest="file", metavar="PATH", type="string", default='assets/image_fbp3d.hv',
                  help="Path to inter_file")

parser.add_option("-v", "--volume_name",
                  action="store", dest="volume_name", type="string", default='sphere22in',
                  help="Path to inter_file")

parser.add_option("-f", "--full",
                  action="store_true", dest="full", default=False,
                  help="print all values")

parser.add_option("-s", "--sigma_gaussian_filter", default=False,
                  action="store", dest="sigma", type="float",
                  help="sigma gaussian filter")

parser.add_option("-m", "--size_median_filter", default=False,
                  action="store", dest="median_size", type="int",
                  help="size median filte")
def interfile_parser(file_name):
    """
    Parse interfile header to dict
    :param file_name: srt interfile path
    :return:
    """
    f = open(file_name, 'r')
    param = {}
    try:
        if file_name.split(".")[-1]== 'raw':
            #parse qetir raw files without header
            param["type"] = np.float32
            param['flip'] = False
            param['scaling_factor_xy'] = 2.5
            param['scaling_factor_z'] = 2.5
            match = re.search("_([1-9]+)x([1-9]+)x([1-9]+)\.raw", file_name)
            param['size'] = (int(match.group(1)), int(match.group(2)), int(match.group(3)))
            param['path_to_data_file'] = file_name
            return param
        else:
            for line in f.readlines():
                matchObj = re.match(r'(.*) := (.*)', line, re.M | re.I)
                if matchObj:
                    param[matchObj.group(1)] = matchObj.group(2)
    except:
        raise Exception("Unsupurted format")
    try:
        param['size'] = (int(param['!matrix size [3]']), int(param['!matrix size [2]']), int(param['!matrix size [1]']))
    except KeyError:
        raise Exception("Bad parsing Matrix size")
    try:
        param['flip'] = bool(param['flip'])
    except:
        param['flip'] = False # add fil for QETRI

    if param['!number format'] == 'float' or param['!number format'] == 'short float':
        if param['!number of bytes per pixel'] == '4':
            param["type"] = np.float32
        else:
            raise Exception("Bad number format")
    else:
        raise Exception("Bad number format")
    try:
        param['scaling_factor_xy'] = float(param['scaling factor (mm/pixel) [1]'])
        param['scaling_factor_z'] = float(param['scaling factor (mm/pixel) [3]'])
    except:
        raise Exception("Bad parsing scaling_factor")
    if param.get('name of data file'):
        param['path_to_data_file'] = path.join(path.dirname(file_name), param['name of data file'])
    elif param.get('!name of data file'):
        param['path_to_data_file'] = path.join(path.dirname(file_name), param['!name of data file'])
    return param


def interfile2array(param):
    """
    conveft interfile to numpy array
    :param param: dict contens interfile poarametrs
    :return:
    """

    f = open(param['path_to_data_file'], 'r')
    v_list = np.fromfile(f, dtype=param['type'])
    f.close()
    resh_arr = np.asarray(v_list).reshape(param['size'])
    if param['flip'] or True:
        resh_arr = resh_arr[:,::-1,::-1] #verfilp for QETIR
    return resh_arr


def get_circle_measure(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z):
    """
    Get mean and std from sphere (x0, y0, z0) [mm]  radius  [mm], scaling_factor mm/pixel, from array with specfic size
    """
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)
    radius_z = int(radius * (1.0 / scaling_factor_z))
    out_list = []
    radius_xy = int(radius * (1 / scaling_factor_xy))

    y, x = np.ogrid[-radius_xy: radius_xy + 1, -radius_xy: radius_xy + 1]
    index = x ** 2 + y ** 2 <= radius_xy ** 2
    out_list.extend(array[z0, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index].flatten())
    # array[z0, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index]=0
    # np.savetxt('array_'+str(z0),array[z0,:,:])
    return np.mean(out_list), np.std(out_list), len(out_list)



def get_sphere_maximum(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z):
    """
    Get max  from sphere (x0, y0, z0) [mm]  radius  [mm], scaling_factor mm/pixel, from array with specfic size
    """
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)
    radius_z = int(radius * (1.0 / scaling_factor_z))
    out_list = []
    z= z0
    if np.abs(np.abs(z0 - z) * scaling_factor_z / radius) < 10e-100:
        radius_xy = radius
    else:
        radius_xy = np.abs(z0 - z) * scaling_factor_z * 1 / tan(asin(np.abs(z0 - z) * scaling_factor_z / radius))
    if isnan(radius_xy):
        radius_xy = radius

    radius_xy = int(radius_xy * (1 / scaling_factor_xy))

    y, x = np.ogrid[-radius_xy: radius_xy + 1, -radius_xy: radius_xy + 1]
    index = x ** 2 + y ** 2 <= radius_xy ** 2
    out_list.extend(array[z, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index].flatten())
    return max(out_list)

def show_xy(_array, r=437.3, title=""):
    """
    Show and save slice with normalize (max value as 1.0)
    """
    pylab.figure()
    ############# code for show image with STIR-TOF ring
    # def cropND(img, bounding):
    #     start = tuple(map(lambda a, da: a//2-da//2, img.shape, bounding))
    #     end = tuple(map(operator.add, start, bounding))
    #     slices = tuple(map(slice, start, end))
    #     return img[slices]

    # # foo[foo == 0] = m
    # size = np.shape(_array)
    # phantom_maximum = max(cropND(_array, (size[0]/2,size[1]/2)).flatten())
    # _array[_array > phantom_maximum]= phantom_maximum

    phantom_maximum = max(_array.flatten())
    _array=(_array-min(_array.flatten()))
    _array = (_array/phantom_maximum)
    pylab.xlabel('X [mm]')
    pylab.ylabel('Y [mm]')
    pylab.xlim([-r, r])
    pylab.imshow(_array, cmap="hot", origin='upper', extent=[-r, r, -r, r])
    fig_name = path.basename(title).split('.')[0]
    pylab.colorbar()
    pylab.savefig(fig_name)
    # pylab.show()

def MSE(A,B):
    return np.square(np.subtract(A, B)).mean()
def norm_image(image):
    phantom_maximum = max(image.flatten())
    image=(image-min(image.flatten()))
    image = (image/phantom_maximum)
    return image

def image_MSE(image_A, image_B):
    image_A_norm = norm_image(image_A)
    image_B_norm = norm_image(image_B)
    return MSE(image_A_norm, image_B_norm)

def raw2array(file_name):
    """
    conveft interfile to numpy array
    :param param: dict contens interfile poarametrs
    :return:
    """

    f = open(file_name, 'r')
    v_list = np.fromfile(f, dtype=np.float32)
    f.close()
    match = re.search("_([1-9]+)x([1-9]+)x([1-9]+)\.raw", file_name)
    size = (int(match.group(1)), int(match.group(2)), int(match.group(3)))
    resh_arr = np.asarray(v_list).reshape(size)
    return resh_arr


def measure(interfile_header, volume_name, full=False, filter_size=None):
    """
    Calculate image quality CRC BV
    """
    parameters = interfile_parser(interfile_header)
    image = interfile2array(parameters) # read intefile
    if filter_size:
        if type(filter_size) == float:
            print("Gaussian filter :%f " % filter_size)
            image = ndimage.gaussian_filter(image, sigma=filter_size)           
        else:
            print("Median filter :%d " % filter_size )
            image = ndimage.median_filter(image, filter_size)

    r = JPET.geometry[volume_name]['setRmax'] * 10 #  JPET has infarmation about NEMA Phantom from GATE

    ##SETS of bacground
    num_backgrounds_sphere = 12

    ## center of bacgraoud ROI like Lech
    background_con = [(111.5105, -1.270829), (72.63656, 56.083368), (7.133568e-15, 81.5), (-67.88933, 59.674647), (-94.70335, 32.849289), (-111.5105, -1.270829), (-116.31280000000001, -39.168228), (-97.61264, -72.413797), (115.48389999999999, -44.667894), (86.48913, -78.478255), (-61.0, -81.5), (-2.5, -81.5)]
    num_backgrounds_sphere = 12

    x, y, z = get_translation(JPET.geometry, volume_name)

    # read image parameters
    scaling_factor_xy = parameters['scaling_factor_xy']
    scaling_factor_z = parameters['scaling_factor_z']
    size_xy = parameters['size'][1]
    size_z = parameters['size'][0]

    back_ground = {'mean': [], "std": []}

    for delta_z in [-20,-10,0,10,20]: # 5 slice -2 cm -1 cm ...
        for i in range(num_backgrounds_sphere):
            m, s, n = get_circle_measure(image, background_con[i][0], background_con[i][1], z+delta_z, r, scaling_factor_xy, scaling_factor_z,
                                    size_xy, size_z)
            # draw_sphere_mm(image, background_con[i][0], background_con[i][1], z+delta_z, r, scaling_factor_xy, scaling_factor_z,
                                    # size_xy, size_z,0)
            back_ground['mean'].append(m)
            back_ground['std'].append(s)

    back_ground_mean = np.mean(back_ground['mean'])
    back_ground_std = np.std(back_ground['mean'])
    k= len(back_ground['mean'])

    # GET ROI

    roi_mean, roi_std, n = get_circle_measure(image, x, y, z, r, scaling_factor_xy, scaling_factor_z, size_xy, size_z)
    # draw_sphere_mm(image, x, y, z, r, scaling_factor_xy, scaling_factor_z, size_xy, size_z,0)



    limit = (scaling_factor_xy*size_xy)/2.0
    slice_z0 = int((z*(1/scaling_factor_z))+(size_z/2)) # 35 mm shpere center in Z
    show_xy(image[slice_z0, :, :], r=limit, title = interfile_header)

    # calculate error
    Delta_bg_mean = back_ground_std/np.sqrt(k)
    Delta_roi_mean = roi_std/np.sqrt(n)
    corr_varince = np.sqrt(2.0/(k-1))
    Delta_SD = back_ground_std * corr_varince
    Delta_CRC = 1/3.0 * np.sqrt( (Delta_roi_mean/back_ground_mean )**2 + (Delta_bg_mean*(roi_mean/(back_ground_mean)**2) )**2)
    Delta_BV = np.sqrt(((back_ground_std*Delta_bg_mean)/(back_ground_mean)**2)**2 + (Delta_SD/back_ground_mean)**2)

    if full: # print std BV and std ROI
        image_ref= raw2array("reference_161x161x161.raw")
        image_ref=norm_image(image_ref)
        image_n= norm_image(image)
        mse = image_MSE(image_n, image_ref)

        bv = (back_ground_std / back_ground_mean)
        bv_error = Delta_BV
        crc = ((roi_mean / back_ground_mean - 1) / 3)
        crc_error = Delta_CRC
        out_str = '{}\t{}\t{}\t{}\t{}\t{}'.format(bv, bv_error, crc,  crc_error, mse, volume_name)
        print(out_str)
    else:
        bv = "%.3f" % (back_ground_std / back_ground_mean)
        bv_error = "%.3f" % Delta_BV
        crc = "%.2f" % ((roi_mean / back_ground_mean - 1) / 3)
        crc_error = "%.2f" % Delta_CRC
        print(' {}            {}'.format('BV', 'CRC'))
        print('{}+-{}  {}+-{}'.format(bv, bv_error, crc, crc_error))





if __name__ == "__main__":
    (options, args) = parser.parse_args()
    file_name = options.file
    volume_name = options.volume_name
    full = options.full #print not rounded values
    sigma = options.sigma
    median_size = options.median_size
    filter_size = sigma if sigma else median_size
    measure(file_name, volume_name, full=full, filter_size=filter_size)
