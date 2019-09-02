from jinja2 import Template
from math import atan,ceil, degrees, radians, tan
# settings
number_of_rings = 100#400
lenght = 50.0  # in cm
inner_ring_diameter = 85.56 # 43.56*2  # in cm
average_depth_of_interaction = 0.95#1.0  # in cm
number_of_detector = 384
FOV = 50.0# in cm

# TOF parameters
tof=True
number_of_tof_bins = 410
size_of_time_bin = 10 # ps
time_resolution = 400 # ps
tof_mashing_factor = 82
#numer of TOF position = number_of_tof_bins / tof_mashing_facto




# outputs file name
sinogram_file_name = "scanner_template.hs"
root_header_file_name = "root_header.hroot"


# calculations

dist = lenght / number_of_rings
non_arc_corrected_bins = int(number_of_detector / 2)
default_bin_size = ((inner_ring_diameter / 2) * 3.14) / number_of_detector

FOV_range = int(int(FOV / lenght * number_of_rings))

size = list(range(number_of_rings - FOV_range+1, number_of_rings + 1))
size.extend(list(range(number_of_rings - 1, number_of_rings - FOV_range, -1)))

diff = list(range(-FOV_range + 1, 1))
diff.extend(list(range(1, FOV_range)))

angle= degrees(atan(FOV/85.56/2.0))
FOV_1 = tan(radians(int(angle)))*(inner_ring_diameter/2.0)
FOV_2 = tan(radians(ceil(angle)))*(inner_ring_diameter/2.0)
FOV1_range = int(int(FOV_1 / lenght * number_of_rings))
FOV2_range = int(int(FOV_2 / lenght * number_of_rings))
# angle= degrees(atan(FOV/85.56/2.0))

print(angle, int(angle), ceil(angle), FOV1_range, FOV2_range, (number_of_rings-FOV1_range+1, number_of_rings-FOV2_range+1))

def list2str(_list):
    """convert list to string fotmat list in STIR"""
    return str(_list).replace(" ", "").replace('[', '{ ').replace(']', ' }')


with open('scanner_template.template') as file_:
    template = Template(file_.read())

with open('root_header.template') as file_:
    template_root = Template(file_.read())


# render templates
output = template.render(number_of_rings=number_of_rings, diff=list2str(diff), size2=list2str(size),
                         dist=dist,
                         size4=len(size), size1=int(number_of_detector / 2),
                         size3=int(number_of_detector / 2),
                         average_depth_of_interaction=average_depth_of_interaction,
                         inner_ring_diameter=inner_ring_diameter,
                         non_arc_corrected_bins=non_arc_corrected_bins,
                         default_bin_size=default_bin_size,
                         number_of_detector=number_of_detector,
                         tof=tof,
                         number_of_tof_bins=number_of_tof_bins,
                         size_of_time_bin=size_of_time_bin,
                         time_resolution=time_resolution,
                         tof_mashing_factor=tof_mashing_factor,
                         )

text_file = open(sinogram_file_name, "w")
text_file.write(output)
text_file.close()
output_root = template_root.render(number_of_rings=number_of_rings, dist=dist,
                                   non_arc_corrected_bins=non_arc_corrected_bins,
                                   default_bin_size=default_bin_size,
                                   number_of_detector=number_of_detector,
                                   inner_ring_diameter=inner_ring_diameter,
                                   average_depth_of_interaction=average_depth_of_interaction,
                                   tof=tof,
                                   number_of_tof_bins=number_of_tof_bins,
                                   size_of_time_bin=size_of_time_bin,
                                   time_resolution=time_resolution,
                                   tof_mashing_factor=tof_mashing_factor,
                                   )

text_file2 = open("root_header.hroot", "w")
text_file2.write(output_root)
text_file2.close()
