#include <osgEarth/MouseCoordsTool>
#include <osgEarth/LatLongFormatter>
#include <osgViewer/Viewer>
#include <osgEarth/MapNode>
#include <osgEarth/EarthManipulator>
#include <osgEarth/NodeUtils>
#include <osgEarth/ExampleResources>
#include <osgEarth/FeatureNode>
#include <osgEarth/TessellateOperator>
#include <osgEarth/SpatialReference>
#include <osgEarth/PlaceNode>
#include <osgDB/ReadFile>
//#include <osg/EllipsoidModel>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/LineWidth>
#include <osg/StateSet>
#include <osg/Material>
#include <osg/Light>
#include <osg/LightSource>
#include <osgViewer/Viewer>
#include <osg/BlendFunc>

#include "mysatellite.h"
#include "coreLib.h"
#include "orbitLib.h"


using namespace osgEarth;
using namespace osgEarth::Util;
using namespace osgEarth::Contrib;
namespace ui = osgEarth::Util::Controls;

        const double IAU_EARTH_ANGULAR_VELOCITY = 7292115.1467e-11; // (rad/sec)

struct App
{
    App() {
        _playing = false;
        _eci = false;
        _show._cone = false;
        _show._cov = false;
        _show._orbit = true;
	_show._label = true;
	_show._scale = 1;
	_show._eleva = 10;
	_show._beamelev = 60;
        _speed = 10;
        readout = new ui::LabelControl();
        readout->setVertAlign(ui::Control::ALIGN_CENTER);
    }

    osg::ref_ptr<PlaceNode> sunPos;
    osg::ref_ptr<PlaceNode> targetPos;
    SkyNode* sky;
    ui::LabelControl* readout;
    ui::LabelControl* lspeed;
    ui::LabelControl* lscale;
    ui::LabelControl* lelev;
    ui::LabelControl* lbeamelev;
    ui::HSliderControl* speed;
    ui::HSliderControl* scale;
    ui::HSliderControl* elev;
    ui::HSliderControl* beamelev;
    ui::CheckBoxControl* eci;
    ui::CheckBoxControl* cone;
    ui::CheckBoxControl* cov;
    ui::CheckBoxControl* orbit;
    ui::CheckBoxControl* lab;
    void play() { _playing = true; }
    void stop() { _playing = false; }

    void tick() {
        if (_playing) {
            TimeStamp t = sky->getDateTime().asTimeStamp() + _speed;
            sky->setDateTime(DateTime(t));
        }
        readout->setText(sky->getDateTime().asRFC1123());
    }
    void setSpeed(){_speed=speed->getValue();
	std::stringstream ss;
	ss<<std::right<<std::setw(4)<<_speed;
	lspeed->setText("speed "+ss.str());}
    void setScale(){_show._scale=scale->getValue();
	std::stringstream ss;
	ss<<std::right<<std::setw(4)<<std::fixed<<std::setprecision(1)<<_show._scale;
	lscale->setText("scale "+ss.str());}
    void setElev(){_show._eleva=elev->getValue();
	std::stringstream ss;
	ss<<std::right<<std::setw(4)<<std::fixed<<std::setprecision(1)<<_show._eleva;
	lelev->setText("ele "+ss.str());}
    void setBeamelev(){_show._beamelev=beamelev->getValue();
	std::stringstream ss;
	ss<<std::right<<std::setw(4)<<std::fixed<<std::setprecision(1)<<_show._beamelev;
	lbeamelev->setText("bele "+ss.str());}
    void setECI(){_eci=eci->getValue();}
    void setCone(){_show._cone=cone->getValue();}
    void setCov(){_show._cov=cov->getValue();}
    void setOrbit(){_show._orbit=orbit->getValue();}
    void setLabel(){_show._label=lab->getValue();}
    bool _eci;
    Showflags _show;

    bool _playing;
    int _speed;
};
OE_UI_HANDLER(setSpeed);
OE_UI_HANDLER(setScale);
OE_UI_HANDLER(setElev);
OE_UI_HANDLER(setBeamelev);
OE_UI_HANDLER(setECI);
OE_UI_HANDLER(setCone);
OE_UI_HANDLER(setCov);
OE_UI_HANDLER(setOrbit);
OE_UI_HANDLER(setLabel);

struct Play : public ui::ControlEventHandler {
    Play(App& app) : _app(app) { }
    void onClick(ui::Control*) { _app.play(); }
    App& _app;
};

struct Stop : public ui::ControlEventHandler {
    Stop(App& app) : _app(app) { }
    void onClick(ui::Control*) { _app.stop(); }
    App& _app;
};

ui::Container* createUI(App& app)
{
    ui::Grid* grid = new ui::Grid();
    grid->setControl(2,0,new ui::ButtonControl("Play", new Play(app)));
    grid->setControl(3,0,new ui::ButtonControl("Stop", new Stop(app)));
    app.lspeed=grid->setControl(4,0,new ui::LabelControl("speed"));
    app.speed=grid->setControl(5,0,new ui::HSliderControl(-100,100,10,new setSpeed(app)));
    app.speed->setWidth(150);
    grid->setControl(6,0,app.readout);
    int col=0;
    app.lscale=grid->setControl(0,col,new ui::LabelControl("model scale"));
    app.scale=grid->setControl(1,col++,new ui::HSliderControl(0.1,5,1,new setScale(app)));
    app.scale->setWidth(100);
    app.lelev=grid->setControl(0,col,new ui::LabelControl("elev thrd"));
    app.elev=grid->setControl(1,col++,new ui::HSliderControl(1,70,10,new setElev(app)));
    app.elev->setWidth(100);
    app.lbeamelev=grid->setControl(0,col,new ui::LabelControl("Beam elev"));
    app.beamelev=grid->setControl(1,col++,new ui::HSliderControl(5,89,60,new setBeamelev(app)));
    app.beamelev->setWidth(100);
    grid->setControl(0,col,new ui::LabelControl("ECI"));
    app.eci=grid->setControl(1,col++,new ui::CheckBoxControl(false,new setECI(app)));
    grid->setControl(0,col,new ui::LabelControl("Beam Cone"));
    app.cone=grid->setControl(1,col++,new ui::CheckBoxControl(false,new setCone(app)));
    grid->setControl(0,col,new ui::LabelControl("Beam Cover"));
    app.cov=grid->setControl(1,col++,new ui::CheckBoxControl(false,new setCov(app)));
    grid->setControl(0,col,new ui::LabelControl("Orbit"));
    app.orbit=grid->setControl(1,col++,new ui::CheckBoxControl(true,new setOrbit(app)));
    grid->setControl(0,col,new ui::LabelControl("Label"));
    app.lab=grid->setControl(1,col++,new ui::CheckBoxControl(true,new setLabel(app)));

    return grid;
}
class CoordinateCollector : public MouseCoordsTool::Callback
{
public:
    GeoPoint currentPoint;
    
    void set(const GeoPoint& coords, osg::View* view, MapNode* mapNode) override
    {
        currentPoint = coords;
    }
    
    void reset(osg::View* view, MapNode* mapNode) override
    {
        currentPoint = GeoPoint(); // 设置为无效点
    }
};

class ClickHandler : public osgGA::GUIEventHandler
{
public:
    ClickHandler(CoordinateCollector* collector,    osg::ref_ptr<PlaceNode> targetPos)
 : _collector(collector), _targetPos(targetPos) {}
    bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
    {
        if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE && 
            ea.getButton() == osgGA::GUIEventAdapter::RIGHT_MOUSE_BUTTON)
        {
            // 从 MouseCoordsTool 获取当前坐标
            const GeoPoint& point = _collector->currentPoint;
            
            if (point.isValid())
            {
                // 打印WGS84坐标
                std::cout << "点击位置坐标: " 
                          << "Lat=" << point.y() << ", "
                          << "Lon=" << point.x() << ", "
                          << "Alt=" << point.z() << "m" << std::endl;
                _targetPos->setPosition(point);
                // 添加标签节点
                //addLabel(point);
            }
            
            return true;
        }
        return false;
    }
    CoordinateCollector* _collector;
    osg::ref_ptr<PlaceNode> _targetPos;

};
int main(int argc, char** argv)
{
    osgEarth::initialize();

    osg::ArgumentParser arguments(&argc,argv);
    vector<cSatellite*> satlist;
    std::string tlefile;
    if(arguments.read("--tle",tlefile))
    {
	readtlefile(tlefile,&satlist);
    }

    osgViewer::Viewer viewer(arguments);
    // Tell the database pager to not modify the unref settings
    viewer.getDatabasePager()->setUnrefImageDataAfterApplyPolicy( false, false );
    // install our default manipulator (do this before calling load)
    osgEarth::Util::EarthManipulator* manip = new EarthManipulator();
    viewer.setCameraManipulator( manip );
//    viewer.setCameraManipulator( new EarthManipulator(arguments) );

    viewer.getCamera()->setClearColor( osg::Vec4(0,0,0,1) );

    App app;

    // load an earth file, and support all or our example command-line options
    // and earth file <external> tags
    osg::Node* node = MapNodeHelper().load( arguments, &viewer );
    if ( node )
    {
        osg::Group* root = new osg::Group();
        root->addChild( node );
        MapNode* mapNode = MapNode::get(node);
	const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
	osgEarth::Viewpoint vp;
	vp.setFocalPoint(GeoPoint(geoSRS,120.0,0.0, ALTMODE_ABSOLUTE));
	vp.range()->set(180000000, Units::METERS);
	manip->setHomeViewpoint( vp );

//osgEarth::Util::SkyNode* skyNode =osgEarth::findTopMostNodeOfType<SkyNode>(node);
//    skyNode->setSimulationTimeTracksDateTime(true);
    osgEarth::Util::SkyOptions skyOptions;
    skyOptions.ambient() = 0.1;
	//skyOptions.coordinateSystem() = osgEarth::Util::SkyOptions::COORDSYS_ECEF;
    skyOptions.quality()=osgEarth::Util::SkyOptions::Quality::QUALITY_DEFAULT;
    osgEarth::Util::SkyNode* skyNode = osgEarth::Util::SkyNode::create(skyOptions);

/*        std::string ext = mapNode->getMapSRS()->isGeographic() ? "sky_simple" : "sky_gl";
        mapNode->addExtension(osgEarth::Extension::create(ext, skyOptions));*/
	skyNode->setLighting(true);
	root->addChild(skyNode);
	app.sky = skyNode;
        const Ephemeris* ephemeris = 0L;
	osg::ref_ptr<osg::Image> mark = osgDB::readRefImageFile("/osgearth/data/placemark32.png");
        if ( app.sky )
        {
            ephemeris = app.sky->getEphemeris();
        app.sunPos = new PlaceNode("Sun", Style(), mark.get());
        app.sunPos->setDynamic(true);
        mapNode->addChild( app.sunPos.get() );

        }

	app.targetPos = new PlaceNode("Target",Style(),mark.get());
        app.targetPos->setDynamic(true);
	app.targetPos->setPosition(GeoPoint(geoSRS, 110, 30));
        mapNode->addChild( app.targetPos.get() );

        ui::ControlCanvas* container = ui::ControlCanvas::getOrCreate(&viewer);
        container->addChild(createUI(app));

	MouseCoordsTool* tool = new MouseCoordsTool(mapNode);
	LatLongFormatter formatter;
	//tool->addCallback(new MouseCoordsLabelCallback(app.target, &formatter));
	viewer.addEventHandler(tool);
	CoordinateCollector* collector = new CoordinateCollector();
	tool->addCallback(collector);

	viewer.addEventHandler(new ClickHandler(collector,app.targetPos));
    // A lat/long SRS for specifying points.
    //const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
    //std::map<int,osg::ref_ptr<osg::Geode> > satobj;
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	std::map<int,SatelliteObj* > satobjmap;
        const DateTime& t = app.sky->getDateTime();
	CelestialBody sun = ephemeris->getSunPosition(t);
	osg::Vec3d sunPos = sun.geocentric;
	sunPos.normalize();
	for(int i=0;i<satlist.size();i++)
    	{
    	    cSatellite* sat = satlist[i];
	    if(sat->Name().find("BEIDOU-3"/*"GPS"*/)!=string::npos){
	    //if(sat->Name().find("POLAR"/*"BEIDOU-3 M1""GPS BIII-6"*/)!=string::npos){
	    //if(sat->Name().find("QZS"/*"GPS BIII-6"*/)!=string::npos){
	    //if(sat->Name().find("GPS BIIF")!=string::npos){
		printf("%s\n",sat->Name().c_str());
		SatelliteObj * satobj = new SatelliteObj(sat,mapNode,geode,t,sunPos,app.targetPos);
		//root->addChild(satobj->boxTransform);
		root->addChild(satobj->pat);
		satobjmap[i]=satobj;
	    }
	}
	//GeodeFinder finder;
    //root->accept(finder);
	root->addChild(geode);


for (const auto& pair : satobjmap)
{
  std::cout << "Key: " << pair.first << ", Value: " << pair.second->sat->Name() << std::endl;
}


        viewer.setSceneData( root );

	//DateTime t0 = app.sky->getDateTime();
        while(!viewer.done())
        {
            viewer.frame();

            if ( ephemeris && app._playing)
            {
                const DateTime& dt = app.sky->getDateTime();
		sun = ephemeris->getSunPosition(dt);
		sunPos = sun.geocentric;
		sunPos.normalize();
		for(const auto & pair : satobjmap){
		    //pair.second->scale = app._scale/2.0;
		    pair.second->showflag = app._show;
//printf("scale %f ele %f beamele %f\n",pair.second->showflag._scale,pair.second->showflag._eleva,pair.second->showflag._beamelev);
		    pair.second->setposition(dt,sunPos);
		}
		if(app._eci){
			//printf("dt %d\n",app._speed);
		    osgEarth::Viewpoint vp0=manip->getViewpoint();
		    auto optionalFocalPoint=vp0.focalPoint();
		    osgEarth::GeoPoint p0=*optionalFocalPoint;
		    p0.x()-=app._speed*OMEGA_E;
		    //vp0.setFocalPoint(GeoPoint(geoSRS,p0.x(),p0.y(), ALTMODE_ABSOLUTE));
		    vp0.setFocalPoint(p0);
		    manip->setViewpoint(vp0);

		    //printf("vp: %lf %lf %lf , %d, %lf\n",p0.x(),p0.y(),p0.z(),app._speed,app._speed*OMEGA_E);
		
		//printf("dt: %d\n",dt.asTimeStamp()-t0.asTimeStamp());
		//t0=dt;

		}
                GeoPoint sunPos1;
                sunPos1.fromWorld(mapNode->getMapSRS(), sun.geocentric);
                sunPos1.alt() = 0.0;
                app.sunPos->setPosition( sunPos1 );
                app.sunPos->setText( "Sun\n");// + llf.format(sunPos1) );

            }

            app.tick();
        }

    }
    return viewer.run();
}
                /*CelestialBody sun = ephemeris->getSunPosition(dt);
                GeoPoint sunPos;
                sunPos.fromWorld(mapNode->getMapSRS(), sun.geocentric);
                sunPos.alt() = 0.0;
                app.sunPos->setPosition( sunPos );
                app.sunPos->setText( "Sun\n" + llf.format(sunPos) );

                CelestialBody moon = ephemeris->getMoonPosition(dt);
                GeoPoint moonPos;
                moonPos.fromWorld(mapNode->getMapSRS(), moon.geocentric);
                moonPos.alt() = 0.0;
                app.moonPos->setPosition( moonPos );
                app.moonPos->setText( "Moon\n" + llf.format(moonPos) );*/
/*class GeodeFinder : public osg::NodeVisitor
{
public:
    GeodeFinder() : osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN) {}

    virtual void apply(osg::Node& node)
    {
        // 如果当前节点是 Geode，处理它
        osg::Geode* geode = dynamic_cast<osg::Geode*>(&node);
        if (geode)
        {
            std::cout << "Found Geode with " << geode->getNumDrawables() << " drawables." << std::endl;

            // 遍历 Geode 中的每个 Drawable
            for (unsigned int i = 0; i < geode->getNumDrawables(); ++i)
            {
                osg::Drawable* drawable = geode->getDrawable(i);
                if (drawable)
                {
                    std::cout << "  Drawable " << i << ": " << drawable->className() << std::endl;
                }
            }
        }

        // 继续遍历子节点
        traverse(node);
    }
};*/
	        //osg::ref_ptr<osg::Geode> geode = createsatellite(sat,mapNode,t);
    /*unsigned int numDrawables = geode->getNumDrawables();
    for (unsigned int i = 0; i < numDrawables; ++i)
    {
        osg::Drawable* drawable = geode->getDrawable(i);
        if (drawable)
        {
            // 处理每个 Drawable 对象
            std::cout << "Drawable " << i << ": " << drawable->className() << std::endl;
            osg::Geometry* geometry = drawable->asGeometry();
            if (geometry)
            {
                // 如果是 Geometry 类型，可以访问其顶点数据等
                osg::Vec3Array* vertices = dynamic_cast<osg::Vec3Array*>(geometry->getVertexArray());
                if (vertices)
                {
                    std::cout << "Geometry has " << vertices->size() << " vertices." << std::endl;
                }
            }
	}
    }*/
    	        //root->addChild(geode);
		//satobj[i]=geode;

    /*//set window title
    osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
    traits->x = 100;
    traits->y = 100;
    traits->width = 800;
    traits->height = 600;
    traits->windowDecoration = true;
    traits->doubleBuffer = true;
    traits->sharedContext = 0;
    traits->windowName = "自定义窗口标题";
    // 创建图形上下文
    osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits);
    if (gc.valid())
    {
        osg::Camera* camera = viewer.getCamera();
        camera->setGraphicsContext(gc);
        camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
        GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
        camera->setDrawBuffer(buffer);
        camera->setReadBuffer(buffer);
    }*/

    /*skyNode->setSunVisible(true);
    skyNode->setMoonVisible(true);
    skyNode->setStarsVisible(true);
    skyNode->setLighting(true);
//    skyNode->setSimulationTimeTracksDateTime(true);*
osg::Light * light = skyNode->getSunLight();
printf("Light: %d\n",light->getLightNum());
osg::Vec4 vec = light->getPosition();
printf("Light pos %lf %lf %lf %lf\n",vec[0],vec[1],vec[2],vec[3]);
vec = light->getAmbient();
printf("Light Ambient %lf %lf %lf %lf\n",vec[0],vec[1],vec[2],vec[3]);
light->setDiffuse( osg::Vec4(.01,.01,.01,0.1) );
vec = light->getDiffuse();
printf("Light Diffuse %lf %lf %lf %lf\n",vec[0],vec[1],vec[2],vec[3]);
osg::Vec3 dir = light->getDirection();
printf("Light Direction %lf %lf %lf\n",dir[0],dir[1],dir[2]);
printf("Light SpotExponent: %f\n",light->getSpotExponent());
printf("Light SpotCutoff: %f\n",light->getSpotCutoff());
//light->setAmbient( osg::Vec4(0.4f, 0.4f, 0.4f ,1.0) );
printf("Light LA %f\n",light->getLinearAttenuation());
light->setLinearAttenuation(0.5);
printf("Light LA %f\n",light->getLinearAttenuation());
printf("Light QA %f\n",light->getQuadraticAttenuation());
light->setQuadraticAttenuation(0.001);
printf("Light QA %f\n",light->getQuadraticAttenuation());

//viewer.setLight(light);
osg::ref_ptr<osg::LightSource> lightSource = new osg::LightSource();
lightSource->setLight(light);
root->addChild(lightSource);*
    osg::ref_ptr<osg::StateSet> stateSet = mapNode->getOrCreateStateSet();
*printf("stateSet %d\n",stateSet->getMode(GL_LIGHTING));
    stateSet->setMode(GL_LIGHTING, osg::StateAttribute::ON);
printf("stateSet GL_LIGHTING %d\n",stateSet->getMode(GL_LIGHTING));
printf("stateSet GL_LIGHT0 %d\n",stateSet->getMode(GL_LIGHT0));
    stateSet->setMode(GL_LIGHT0, osg::StateAttribute::ON);
printf("stateSet GL_LIGHT0 %d\n",stateSet->getMode(GL_LIGHT0));*
osg::ref_ptr<osg::Material> material = new osg::Material();
material->setAmbient(osg::Material::FRONT, osg::Vec4(0.1, 0.1, 0.1, 1.0)); 
material->setDiffuse(osg::Material::FRONT, osg::Vec4(0.15f, 0.15f, 0.15f, 1.0f));
    stateSet->setAttributeAndModes(material, osg::StateAttribute::ON);

//skyNode->getDefaultLight()->setAmbient(osg::Vec4(0.2f, 0.2f, 0.2f, 1.0f));
//skyNode->setTimeControlEnabled(false);*/

/*
void readtlefile(string tlefile, vector<cSatellite*>* satlist){
        std::ifstream fin(tlefile.c_str());
        std::string linename, line1, line2;
        int noradid;
        while(!fin.eof())
        {
            std::string tmp;
            std::getline(fin,tmp);
            int flag=osgEarth::as<int>(tmp.substr(0,1),99);
            //std::cout<<flag<<std::endl;
            if(!tmp.empty())
                tmp.pop_back();
            if(flag==1)
                line1=tmp;
            else if(flag==2){
                line2=tmp;
                noradid=osgEarth::as<int>(tmp.substr(1,6),0);
                //printf("noradid %d\n",noradid);
                if(linename.empty())
                {
                /*    std::cout<<linename;
                    linename="";
                }else{*
                    linename=line1.substr(1,7);
                }
                std::cout<<noradid<<std::endl;
                cTle tle(linename, line1, line2);
                cSatellite* sat=new cSatellite(tle);
                time_t ep = sat->Orbit().Epoch().ToTime();
                printf("sat %s ep %s %ld\n",sat->Name().c_str(),ctime(&ep),ep);
                satlist->push_back(sat);
            }else
                linename=tmp;
        }
    printf("read %d satellites\n",satlist->size());
} 
osg::ref_ptr<osg::Geometry> createColoredBox(osg::Vec3 center, osg::Vec3 scale)
{
    // 创建顶点数组
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    // 定义盒子的8个顶点
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(), -0.5f*scale.y()+center.y(), -0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(), -0.5f*scale.y()+center.y(), -0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(),  0.5f*scale.y()+center.y(), -0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(),  0.5f*scale.y()+center.y(), -0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(), -0.5f*scale.y()+center.y(),  0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(), -0.5f*scale.y()+center.y(),  0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(),  0.5f*scale.y()+center.y(),  0.5f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(),  0.5f*scale.y()+center.y(),  0.5f*scale.z()+center.z()));

    vertices->push_back(osg::Vec3( 0.0f*scale.x()+center.x(),  0.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(),  1.0f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(),  2.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(),  2.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(),  1.0f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));

    vertices->push_back(osg::Vec3( 0.0f*scale.x()+center.x(), -0.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(), -1.0f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3( 0.5f*scale.x()+center.x(), -2.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(), -2.5f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));
    vertices->push_back(osg::Vec3(-0.5f*scale.x()+center.x(), -1.0f*scale.y()+center.y(),  0.0f*scale.z()+center.z()));

    // 创建颜色数组
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    // 为每个面定义不同的颜色
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f)); // 前面 - 红色
    colors->push_back(osg::Vec4(0.0f, 1.0f, 0.0f, 1.0f)); // 后面 - 绿色
    colors->push_back(osg::Vec4(0.0f, 0.0f, 1.0f, 1.0f)); // 顶面 - 蓝色
    colors->push_back(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f)); // 底面 - 黄色
    colors->push_back(osg::Vec4(1.0f, 0.0f, 1.0f, 1.0f)); // 左面 - 紫色
    colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f)); // 右面 - 青色

    colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f)); // 右面 - 青色
    colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f)); // 右面 - 青色

    // 创建法线数组
    osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array;
    normals->push_back(osg::Vec3(0.0f, 0.0f, 1.0f));
    normals->push_back(osg::Vec3(0.0f, 0.0f, -1.0f));
    normals->push_back(osg::Vec3(0.0f, 1.0f, 0.0f));
    normals->push_back(osg::Vec3(0.0f, -1.0f, 0.0f));
    normals->push_back(osg::Vec3(-1.0f, 0.0f, 0.0f));
    normals->push_back(osg::Vec3(1.0f, 0.0f, 0.0f));

    // 创建 Geometry 对象
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;

    // 设置顶点数据
    geometry->setVertexArray(vertices);

    // 设置颜色数据
    geometry->setColorArray(colors);
    geometry->setColorBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);

    // 设置法线数据
    geometry->setNormalArray(normals);
    geometry->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);

    // 定义六个面的顶点索引
    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS, 4, 4)); // 前面
    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS, 0, 4)); // 后面

    osg::ref_ptr<osg::DrawElementsUInt> drawElements = new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS);
    osg::ref_ptr<osg::DrawElementsUInt> drawElements1 = new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS);
    osg::ref_ptr<osg::DrawElementsUInt> drawElements2 = new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS);
    osg::ref_ptr<osg::DrawElementsUInt> drawElements3 = new osg::DrawElementsUInt(osg::PrimitiveSet::QUADS);
    osg::ref_ptr<osg::DrawElementsUInt> drawElements4 = new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON);
    osg::ref_ptr<osg::DrawElementsUInt> drawElements5 = new osg::DrawElementsUInt(osg::PrimitiveSet::POLYGON);
drawElements->push_back(2);
drawElements->push_back(3);
drawElements->push_back(7);
drawElements->push_back(6);
drawElements1->push_back(0);
drawElements1->push_back(1);
drawElements1->push_back(5);
drawElements1->push_back(4);
drawElements2->push_back(0);
drawElements2->push_back(3);
drawElements2->push_back(7);
drawElements2->push_back(4);
drawElements3->push_back(1);
drawElements3->push_back(2);
drawElements3->push_back(6);
drawElements3->push_back(5);

drawElements4->push_back(8);
drawElements4->push_back(9);
drawElements4->push_back(10);
drawElements4->push_back(11);
drawElements4->push_back(12);

drawElements5->push_back(13);
drawElements5->push_back(14);
drawElements5->push_back(15);
drawElements5->push_back(16);
drawElements5->push_back(17);


geometry->addPrimitiveSet(drawElements);
geometry->addPrimitiveSet(drawElements1);
geometry->addPrimitiveSet(drawElements2);
geometry->addPrimitiveSet(drawElements3);
geometry->addPrimitiveSet(drawElements4);
geometry->addPrimitiveSet(drawElements5);

    return geometry;
}
vector<double> SatEcf(const cSatellite& sat, time_t t)
{
    cJulian cjt(t);
    cEciTime sateci = sat.PositionEci(cjt);
    cGeo satgeo(sateci, cjt);
    //printf("%s", ctime(&t));
    //printf("%d %f %f %f\n", t, satgeo.LongitudeDeg(), satgeo.LatitudeDeg(), satgeo.AltitudeKm());
    vector<double>llh;
    llh.push_back(satgeo.LongitudeDeg());
    llh.push_back(satgeo.LatitudeDeg());
    llh.push_back(satgeo.AltitudeKm());
    return llh;
}
*/
