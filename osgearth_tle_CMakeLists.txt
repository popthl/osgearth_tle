INCLUDE_DIRECTORIES(${OSG_INCLUDE_DIRS} )
SET(TARGET_LIBRARIES_VARS OSG_LIBRARY OSGDB_LIBRARY OSGUTIL_LIBRARY OSGVIEWER_LIBRARY OPENTHREADS_LIBRARY)

SET(TARGET_SRC osgearth_tle.cpp mysatellite.cpp)

#### end var setup  ###
SETUP_APPLICATION(osgearth_tle)
set(INCLUDEDIR ${PROJECT_SOURCE_DIR}/orbittools/core ${PROJECT_SOURCE_DIR}/orbittools/orbit)
MESSAGE("mytest includedir:" ${INCLUDEDIR})
INCLUDE_DIRECTORIES(${INCLUDEDIR})
#LINK_EXTERNAL(orbitTools)
MESSAGE(mytest_osgearth_tle)
#target_link_options(mytest_osgearth_mytest PRIVATE -L/usr/lib -lorbitTools)
target_link_libraries(mytest_osgearth_tle /usr/lib/liborbitTools.so)
