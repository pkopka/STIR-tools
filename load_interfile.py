import numpy as np
import pylab
import re
from os import path
from math import asin, tan, isnan
from image_tools import draw_sphere_mm, get_translation
from gate_macro_parser.phantom_parser import JPET

from optparse import OptionParser


parser = OptionParser()

parser.add_option("-i", "--inter_file",
                  action="store", dest="file", metavar="PATH", type="string", default='assets/image_fbp3d.hv',
                  help="Path to inter_file")

parser.add_option("-v", "--volume_name",
                  action="store", dest="volume_name", type="string", default='sphere22in',
                  help="Path to inter_file")



def interfile_parser(file_name):
    """
    Parse interfile header to dict
    :param file_name: srt interfile path
    :return:
    """
    f = open(file_name, 'r')
    param = {}
    for line in f.readlines():
        matchObj = re.match(r'(.*) := (.*)', line, re.M | re.I)
        if matchObj:
            param[matchObj.group(1)] = matchObj.group(2)
    try:
        param['size'] = (int(param['!matrix size [3]']), int(param['!matrix size [2]']), int(param['!matrix size [1]']))
    except KeyError:
        raise Exception("Bad parsing Matrix size")

    if param['!number format'] == 'float':
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
    param['path_to_data_file'] = path.join(path.dirname(file_name), param['name of data file'])
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
    return resh_arr


def get_sphere_measure(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z):
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)

    radius_z = int(radius * (1.0 / scaling_factor_z))
    out_list = []
    for z in range(z0 - radius_z, z0 + radius_z + 1):
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
    return np.mean(out_list), np.std(out_list)


def show_xy(_array, r=437.3):
    pylab.figure()

    # normalize
    # mini  =min(_array.flatten())
    # _array=(_array-min(_array.flatten()))
    # _array = _array/max(_array.flatten())*4
    pylab.xlabel('X [mm]')
    pylab.ylabel('Y [mm]')
    pylab.xlim([-r, r])

    pylab.imshow(_array, cmap="hot", origin='upper', extent=[-r, r, -r, r])
    pylab.colorbar()

    pylab.show()


def measure(interfile_header, volume_name):
    parameters = interfile_parser(interfile_header)
    image = interfile2array(parameters)
    r = JPET.geometry[volume_name]['setRmax'] * 10

    ##SETS of bacground
    num_backgrounds_sphere = 12
    R = 85  # distans from center phantom to center background sphere cm
    z = 37.0  #cm

    scaling_factor_xy = parameters['scaling_factor_xy']
    scaling_factor_z = parameters['scaling_factor_z']
    size_xy = parameters['size'][1]
    size_z = parameters['size'][0]

    alfa_n = np.linspace(0.0, 2 * np.pi, num=num_backgrounds_sphere)
    back_ground = {'mean': [], "std": []}
    for alfa in alfa_n:
        m, s = get_sphere_measure(image, R * np.cos(alfa), R * np.sin(alfa), z, r, scaling_factor_xy, scaling_factor_z,
                                  size_xy, size_z)
        back_ground['mean'].append(m)
        back_ground['std'].append(s)
        # check ROI background sphere
        # draw_sphere_mm(image, R * np.cos(alfa), R * np.sin(alfa), z, r, scaling_factor_xy, scaling_factor_z, size_xy,
        #                size_z, 0)

    back_ground_mean = np.mean(back_ground['mean'])
    back_ground_std = np.mean(back_ground['std'])

    x, y, z = get_translation(JPET.geometry, volume_name)
    roi_mean, roi_std = get_sphere_measure(image, x, y, z, r, scaling_factor_xy, scaling_factor_z, size_xy, size_z)


    # check ROI
    # draw_sphere_mm(image, x, y, z, r, scaling_factor_xy, scaling_factor_z, size_xy, size_z, 0)
    # print(min(image.flatten()), 'min')
    # show_xy(image[57, :, :], r=447.76)
    # pylab.show()
    # image.astype(np.float32).tofile('output/test.v')

    bv = "%.2f" % (back_ground_std / back_ground_mean)
    crc = "%.2f" % ((roi_mean / back_ground_mean - 1) / 3)
    print('{}  {}'.format('BV', 'CRC'))
    print('{}  {}'.format(bv, crc))





if __name__ == "__main__":
    (options, args) = parser.parse_args()
    file_name = options.file
    volume_name = options.volume_name
    measure(file_name, volume_name)
