#include <osgEarth/TessellateOperator>
#include <osgEarth/SpatialReference>
#include <osgEarth/MapNode>
#include <osgEarth/LabelNode>
#include <osgEarth/PlaceNode>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/LineWidth>
#include <osg/StateSet>
#include <osg/Material>
#include <osg/PositionAttitudeTransform>
#include <osgEarth/FeatureNode>

#include "coreLib.h"
#include "orbitLib.h"

#define OMEGA_E 0.00415528985974925 //  deg/s   1.00273790934; // earth rotation per sideral day
#define RADIUS  500000.0f

using namespace osgEarth;
using namespace osgEarth::Util;

struct Showflags{
    bool _cone;
    bool _cov;
    bool _orbit;
    bool _label;
    double _scale;
    double _eleva;
    double _beamelev;
};
class SatelliteObj{
public:
    SatelliteObj(cSatellite* sattle,MapNode* mapnode,osg::ref_ptr<osg::Geode> geode, DateTime t0, osg::Vec3d sunPos,osg::ref_ptr<PlaceNode> targetPos);
    int setposition(DateTime t, osg::Vec3d sunPos);
    void setatt(osg::ref_ptr<osg::Vec3Array> vertices, osg::Vec3d sunPos);
    void setvisible(bool vis);
    cSatellite* sat;
    osg::ref_ptr<osg::Geode> satgroup;
    osg::ref_ptr<PlaceNode> target;
    MapNode* mapNode;
    osg::ref_ptr<osg::Geometry> cov;//sat beam cone
    osg::ref_ptr<osg::Geometry> orbit;
    osgEarth::FeatureNode* featureNode;//sat beam coverage
    osg::ref_ptr<osg::Node> box;
    osg::ref_ptr<osg::MatrixTransform> wing;
    LabelNode* label;
    osg::ref_ptr<osg::PositionAttitudeTransform> pat;
    Showflags showflag;
    bool _vis;
};

void readtlefile(string tlefile, vector<cSatellite*>* satlist);
