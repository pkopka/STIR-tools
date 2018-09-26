
STIR - tools
============
Useful tools for STIR users to loading and generating Interfiles.



image tools
-----------

image_tools generation the attenumation map as interfile. Creating an image requires giving the size in pixels and scaling_factor in mm/pixel.


```
python image_tools -out attenumation_map.v --size_xy 512 --size_z 99 --scaling_factor_xy 1.74 --scaling_factor_z 5.0
```



load interfile
--------------
laod_interfile reads interfile header file and returns statistical parameters for the selected volume described Gate macro with phantom geometry.

```
python load_interfile --inter_file image.hv --volume_name sphere22in
```


Reqrements
-----------
* python 3.4