from jinja2 import Template
from math import atan,ceil, degrees, radians, tan
# settings
number_of_rings = 500
lenght = 200.0  # in cm
inner_ring_diameter = 85.56   # in cm
average_depth_of_interaction = 0.95  # in cm
number_of_detector = 384
FOV = 200# in cm


# outputs file name
sinogram_file_name = "nema_quality.hs"
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




def list2str(_list):
    """convert list to string fotmat list in STIR"""
    return str(_list).replace(" ", "").replace('[', '{ ').replace(']', ' }')


with open('nema.template') as file_:
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
                         number_of_detector=number_of_detector
                         )

text_file = open(sinogram_file_name, "w")
text_file.write(output)
text_file.close()
output_root = template_root.render(number_of_rings=number_of_rings, dist=dist,
                                   non_arc_corrected_bins=non_arc_corrected_bins,
                                   default_bin_size=default_bin_size,
                                   number_of_detector=number_of_detector,
                                   inner_ring_diameter=inner_ring_diameter,
                                   average_depth_of_interaction=average_depth_of_interaction, )

text_file2 = open("root_header.hroot", "w")
text_file2.write(output_root)
text_file2.close()
