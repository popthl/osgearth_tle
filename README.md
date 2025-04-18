# osgearth_tle
draw satellites in osgearth

use osgearth3.2 docker image
use Michael F. Henry's Zeptomoby/OrbitTools to dealwith TLE satellites

in wsl run docker:
docker run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v /dev/dri:/dev/dri -v /home/thl/osgearth_cache:/osgearth_cache -e LIBGL_ALWAYS_SOFTWARE=1 -e MESA_GL_VERSION_OVERRIDE=3.3 -e GDAL_DATA=/usr/share/gdal --privileged -v /home/thl/osgearthproj:/osgearth/proj pelicanmapping/osgearth /bin/bash

in the image
cd /osgearth
mkdir proj
edit /osgearth/CMakeList.txt
ADD_SUBDIRECTORY(proj)

cmake .
make
bin/osgearth_tle tests/simple.earth --tle /osgearth_cache/tle.txt

![image](https://github.com/user-attachments/assets/1567eab8-9062-45b2-997d-d3249bc549f2)

