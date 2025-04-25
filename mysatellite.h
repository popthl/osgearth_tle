#include <osgEarth/TessellateOperator>
#include <osgEarth/SpatialReference>
#include <osgEarth/MapNode>
#include <osgEarth/LabelNode>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/LineWidth>
#include <osg/StateSet>
#include <osg/Material>
#include <osg/PositionAttitudeTransform>

#include "coreLib.h"
#include "orbitLib.h"

#define OMEGA_E 0.00415528985974925 //  deg/s   1.00273790934; // earth rotation per sideral day
#define RADIUS  500000.0f

using namespace osgEarth;
using namespace osgEarth::Util;

class SatelliteObj{
public:
    SatelliteObj(cSatellite* sattle,MapNode* mapnode,osg::ref_ptr<osg::Geode> geode, time_t t0, osg::Vec3d sunPos);
    void setposition(time_t t, osg::Vec3d sunPos);
    void setatt(osg::ref_ptr<osg::Vec3Array> vertices, osg::Vec3d sunPos);
    cSatellite* sat;
    osg::ref_ptr<osg::Geode> satgroup;
    MapNode* mapNode;
    osg::ref_ptr<osg::Geometry> orbit;
    osg::ref_ptr<osg::Node> box;
    osg::ref_ptr<osg::MatrixTransform> wing;
    LabelNode* label;
    osg::ref_ptr<osg::PositionAttitudeTransform> pat;
    double scale;
};

osg::ref_ptr<osg::Geometry> createColoredBox(osg::Vec3 center, osg::Vec3 scale);
void readtlefile(string tlefile, vector<cSatellite*>* satlist);
vector<double> SatEcf(const cSatellite& sat, time_t t);
osg::ref_ptr<osg::Geode> createsatellite(cSatellite* sat,osgEarth::MapNode* mapNode,time_t t0);
