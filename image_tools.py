from math import asin, tan, isnan
from optparse import OptionParser

import numpy as np

parser = OptionParser()

parser.add_option("-o", "--out",
                  action="store", dest="file", metavar="PATH", type="string", default='output/att_map_test.v',
                  help="Path out file")
parser.add_option("-x", "--size_xy",
                  action="store", dest="size_xy", type="int", default=512,
                  help="Image size for X and Y dimension ")

parser.add_option("-z", "--size_z",
                  action="store", dest="size_z", type="int", default=99,
                  help="Image size for Z dimension ")

parser.add_option("-s", "--scaling_factor_xy",
                  action="store", dest="scaling_factor_xy", type="float", default=1.74908,
                  help="Scaling Factor for X Y  dimension in mm/pixel")

parser.add_option("-f", "--scaling_factor_z",
                  action="store", dest="scaling_factor_z", type="float", default=5.0,
                  help="Scaling Factor for X Y  dimension in mm/pixel")


def draw_cylinder_mm(array, x0, y0, z0, radius, height, scaling_factor_xy, scaling_factor_z, size_xy, size_z, mu=0,
                     half=False):
    """
    Draw cylinder in numpy array.

    :param array: image numpy array
    :param x0:
    :param y0: coordinates of the center in mm
    :param z0:
    :param radius:
    :param height: cilinder parameters in mm
    :param scaling_factor_xy:
    :param scaling_factor_z: scaling_factor mm/pixel
    :param size_xy: *array* size X
    :param size_z: *array* size Z
    :param mu: attenuation coefficient in cm-1
    :param half: *bool* True if draw semicircle
    :return:
    """

    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)
    radius = int(radius * (1 / scaling_factor_xy))
    h = int((height * 1 / scaling_factor_z) / 2)  # h/2

    y, x = np.ogrid[-radius: radius + 1, -radius: radius + 1]
    if half == True:

        y, x = np.ogrid[-radius:0 + 1, -radius: radius + 1]  # draw bottom half cylinder
        index = x ** 2 + y ** 2 <= radius ** 2
        falses = np.zeros((radius, 2 * radius + 1), dtype=bool)
        index = np.concatenate((index, falses), axis=0)
    else:
        index = x ** 2 + y ** 2 <= radius ** 2
    for x in range(z0 - h, z0 + h + 1):
        array[x, y0 - radius:y0 + radius + 1, x0 - radius:x0 + radius + 1][index] = mu


def draw_sphere_mm(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z, mu=0):
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)

    radius_z = int(radius * (1.0 / scaling_factor_z))

    for z in range(z0 - radius_z, z0 + radius_z + 1):
        if np.abs(np.abs(z0 - z) * scaling_factor_z / radius)< 10e-100:
            radius_xy = radius
        else:
            radius_xy = np.abs(z0 - z) * scaling_factor_z * 1 / tan(asin(np.abs(z0 - z) * scaling_factor_z / radius))
        if isnan(radius_xy):
            radius_xy = radius

        radius_xy = int(radius_xy * (1 / scaling_factor_xy))

        y, x = np.ogrid[-radius_xy: radius_xy + 1, -radius_xy: radius_xy + 1]
        index = x ** 2 + y ** 2 <= radius_xy ** 2
        array[z, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index] = mu


def get_sphere_mm(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z, mu=0):
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)
    radius_z = int(radius * (1.0 / scaling_factor_z))
    out_list = []
    for z in range(z0 - radius_z, z0 + radius_z + 1):
        radius_xy = np.abs(z0 - z) * scaling_factor_z * 1 / tan(asin(np.abs(z0 - z) * scaling_factor_z / radius))
        if isnan(radius_xy):
            radius_xy = radius
        radius_xy = int(radius_xy * (1 / scaling_factor_xy))

        y, x = np.ogrid[-radius_xy: radius_xy + 1, -radius_xy: radius_xy + 1]
        index = x ** 2 + y ** 2 <= radius_xy ** 2
        out_list.extend(array[z, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index].flatten())
    return np.mean(out_list), np.std(out_list)


def draw_circle_mm(array, x0, y0, z0, radius, scaling_factor_xy, scaling_factor_z, size_xy, size_z, mu=0):
    def cm2pix(dm, dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / dx)))

    x0 = cm2pix(x0)
    y0 = cm2pix(-y0)
    z0 = cm2pix(z0, scaling_factor_z, size_z)
    radius_z = int(radius * (1.0 / scaling_factor_z))

    for z in range(z0 - radius_z, z0 + radius_z + 1):
        radius_xy = np.abs(z0 - z) * scaling_factor_z * 1 / tan(asin(np.abs(z0 - z) * scaling_factor_z / radius))
        if isnan(radius_xy):
            radius_xy = radius

        radius_xy = int(radius_xy * (1 / scaling_factor_xy))

        y, x = np.ogrid[-radius_xy: radius_xy + 1, -radius_xy: radius_xy + 1]
        index = x ** 2 + y ** 2 == radius_xy ** 2
        array[z, y0 - radius_xy:y0 + radius_xy + 1, x0 - radius_xy:x0 + radius_xy + 1][index] = mu


def draw_box_mm(array, x0, y0, z0, dx, dy, dz, scaling_factor_xy, scaling_factor_z, size_xy, size_z, mu=0):
    def cm2pix(dm, _dx=scaling_factor_xy, size=size_xy):
        return int(size / 2 + (dm * (1 / _dx)))

    draw_box(array, cm2pix(x0), cm2pix(-y0), cm2pix(z0, scaling_factor_z, size_z), dx * 1 / scaling_factor_xy,
             dy * 1 / scaling_factor_xy, dz * 1 / scaling_factor_z,
             mu)  # -y ishow images


def draw_box(array, x0, y0, z0, dx, dy, dz, mu):
    array[int((z0 - dz / 2) + 1):int((z0 + dz / 2) + 1), int(y0 - dy / 2):int(y0 + dy / 2 + 1),
    int(x0 - dx / 2):int(x0 + dx / 2 + 1)] = mu


def get_translation(_dict, volume):
    """Change units cm -> mm, roate image (imshow)"""
    return -10 * _dict[volume]['setTranslation'][0], 10 * _dict[volume]['setTranslation'][1], 10 * \
           _dict[volume]['setTranslation'][2]

def draw_geometry(array, geometry, scaling_factor_xy, scaling_factor_z):
    """

    :param array: image numpy array
    :param geometry: python model geometry *dict*
    :param scaling_factor_xy:
    :param scaling_factor_z: scaling_factor mm/pixel
    :return:
    """



    def set_material_val(material):
        """setattenuation coefficient in cm-1"""
        if material == 'Water':
            return 0.096
        elif material == 'Air':
            return 0
        elif material == 'Plastic':
            return 0.097
        elif material == 'Abstract':
            return 0
        else:
            print(material)
            raise Exception('Not ivalid material')

    order_keys = sorted(geometry.keys(), key=lambda x: geometry[x]['factor'])  # draw parents value first
    size_z, size, size_y = np.shape(array)
    for volume in order_keys:
        mat_val = set_material_val(geometry[volume]['material'])
        x, y, z = get_translation(geometry, volume)

        if geometry[volume]['type'] == 'box':
            dx = 10 * geometry[volume]['setXLength']
            dy = 10 * geometry[volume]['setYLength']
            dz = 10 * geometry[volume]['setZLength']
            # print(volume, (x,y,z), dx, dy, dz  ,geometry[volume]['type'])
            draw_box_mm(array, x, y, z, dx, dy, dz, scaling_factor_xy, scaling_factor_z, size, size_z, mat_val)
        elif geometry[volume]['type'] == 'sphare':

            r = 10 * geometry[volume]['setRmax']
            # print(volume, (x,y,z), r, geometry[volume]['type'])
            draw_sphere_mm(array, x, y, z, r, scaling_factor_xy, scaling_factor_z, size, size_z, mat_val)
        elif geometry[volume]['type'] == 'cylinder':
            r = 10 * geometry[volume]['setRmax']
            h = 10 * geometry[volume]['setHeight']
            if geometry[volume].get('setDeltaPhi'):
                half = geometry[volume]['setDeltaPhi'] == 180.0 and geometry[volume]['setPhiStart'] == 0.0
            else:
                half = False
            # print(volume, (x,y,z), r, h ,geometry[volume]['type'], half)
            draw_cylinder_mm(array, x, y, z, r, h, scaling_factor_xy, scaling_factor_z, size, size_z, mat_val,
                             half=half)
        else:
            # print(geometry[volume]['type'])
            raise Exception("Not invalid shape type")


if __name__ == "__main__":
    (options, args) = parser.parse_args()

    print("Import Phantom parser")
    from gate_macro_parser.phantom_parser import JPET
    print("Load geometry to python model...")


    scaling_factor_xy = options.scaling_factor_xy
    scaling_factor_z = options.scaling_factor_z
    size_xy = options.size_xy
    size_z = options.size_z
    array = np.zeros((size_z, size_xy, size_xy))

    print("Create phantom...")
    draw_geometry(array, JPET.geometry, scaling_factor_xy, scaling_factor_z)

    filename = options.file
    print("Save image to file %s ..." % filename)
    array.astype(np.float32).tofile(filename)
    print("Image save successfully")
