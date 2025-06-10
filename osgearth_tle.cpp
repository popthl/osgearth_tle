//#include <osgEarthDrivers/sky_simple/SimpleSkyNode>
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
#include <regex>

using namespace osgEarth;
using namespace osgEarth::Util;
using namespace osgEarth::Contrib;
namespace ui = osgEarth::Util::Controls;

const double IAU_EARTH_ANGULAR_VELOCITY = 7292115.1467e-11; // (rad/sec)

// 前向声明
class TreeNode;

// 点击事件处理器
class TreeNodeClickHandler : public ControlEventHandler {
public:
    TreeNodeClickHandler(TreeNode* node) : _node(node) {}
    virtual void onClick(Control* control) override;
    
private:
    osg::observer_ptr<TreeNode> _node;
};


// 树节点类（必须继承自容器控件）
class TreeNode : public Grid { // Grid 是容器控件，支持 addControl()
public:
    TreeNode(const std::string& name, int level= -1, SatelliteObj * satobj = nullptr, bool expanded = false, TreeNode* parent = nullptr);
    
    void addChild(TreeNode* child);
    void updateparent(TreeNode* parent);
    void toggle();
    
    bool isChecked() const { return _checkbox->getValue(); }
    void setChecked(bool checked) { _checkbox->setValue(checked); }
    void setsatidx(int idx){_satidx=idx;}
    int getsatidx(){return _satidx;}
    SatelliteObj* getsatobj(){return _satobj;}
    bool getexpanded(){return _expanded;}
    // 设置展开状态
    void setExpanded(bool expanded) {
        _expanded = expanded;
        _toggle->setText(expanded ? "-" : "+");
        _children->setVisible(expanded);
    }
    // 设置节点名称
    void setName(const std::string& name) {
        if (_label.valid()) _label->setText(name);
    }
    
    // 获取节点名称
    std::string getName() const {
        return _label.valid() ? _label->text() : "";
    }
    VBox* getChildren(){return _children;}
    TreeNode* getparent(){return _parent;}
private:
    std::string _name;
    bool _expanded;
    osg::ref_ptr<ButtonControl> _toggle;
    osg::ref_ptr<LabelControl> _label;
    osg::ref_ptr<VBox> _children;
    TreeNode* _parent;
    CheckBoxControl* _checkbox;
    int _childnum;
    int _satidx;
    SatelliteObj * _satobj;
    /*ButtonControl* _toggle;
    LabelControl* _label;
    VBox* _children;*/ // 子容器必须是 VBox/HBox/Grid 等容器类型
};

//--------------------------------------------------
// 实现部分
//--------------------------------------------------
void TreeNodeClickHandler::onClick(Control* control) {
    //std::cout<<_node->getName()<<"onClick"<<std::endl;
    if (_node.valid()) {
        _node->toggle();
    }
}
// 自定义复选框事件处理器
class CheckboxHandler : public ControlEventHandler {
public:
    CheckboxHandler(TreeNode* node) : _node(node) {}
    
    virtual void onValueChanged(Control* control, bool value) override {
        if (_node.valid()) {
            // 示例：打印节点选中状态
            TreeNode* parent = _node->getparent();
            /*if(parent!=nullptr){
                OE_NOTICE<<"parent"<<parent->getName();
            }
            OE_NOTICE << "Node [" << _node->getName() 
                      << "] checked: " << (value ? "true" : "false") 
                      <<_node->getsatidx()<< std::endl;*/
            if(_node->getsatidx()>=0){
                //OE_NOTICE<<"sat "<<_node->getsatidx()<< std::endl;
                _node->getsatobj()->setvisible(value);
            }else{
                VBox* children = _node->getChildren();
                if(children->getNumChildren()>1){
                    //OE_NOTICE << "children num"<<children->getNumChildren()<<std::endl;
                    for(int i=0;i<children->getNumChildren();i++){
                        TreeNode* child = dynamic_cast<TreeNode*>(children->getChild(i));
                        if(child!=nullptr){
                            //OE_NOTICE<<"child node "<<child->getName()<<child->getsatidx()<<std::endl;
                            child->setChecked(value);
                        }
                    }
                }
            }
        }
    }
    
private:
    osg::observer_ptr<TreeNode> _node;
};
TreeNode::TreeNode(const std::string& name, int level, SatelliteObj * satobj, bool expanded, TreeNode* parent)
    : _name(name), _satidx(level), _satobj(satobj),  _expanded(expanded), _parent(parent)
{
    _expanded = true;
    // 创建展开按钮
    //_toggle = new ButtonControl(_expanded ? "-" : "+");
    _toggle = new ButtonControl("o");
    _toggle->setPadding(5);
    _toggle->addEventHandler(new TreeNodeClickHandler(this));
    
    // 创建标签
    _label = new LabelControl(_name);
    _label->setPadding(5);
    
    // 配置网格布局（2列：按钮+标签）
    setPadding(2);
    setChildSpacing(5);
    _checkbox = new CheckBoxControl();
    _checkbox->addEventHandler(new CheckboxHandler(this)); // ✅ 绑定处理器
    HBox* hb = new HBox();
    hb->addControl(_checkbox);
    hb->addControl(_label.get());
    setControl(1,0,hb);
    setControl(0,0,_toggle.get());  // ✅ Grid 是容器，支持 addControl()
    
    // 创建子节点容器（必须是容器类型）
    _children = new VBox();
    _children->setVisible(_expanded);
    setControl(1,1,_children.get()); // ✅ 正确添加到 Grid 容器
    _childnum = 0;
}

void TreeNode::updateparent(TreeNode* parent) {
    if(parent->getparent()!=nullptr){
        updateparent(parent->getparent());
    }
    parent->_childnum++;
    parent->setName(parent->_name+"("+std::to_string(parent->_childnum)+")");
    
}
void TreeNode::addChild(TreeNode* child) {
    child->_parent = this;
    if(child->getsatidx()>=0){
        _childnum++;
        if(_parent!=nullptr){
            updateparent(_parent);
        }
    }
    setName(_name+"("+std::to_string(_childnum)+")");
    _children->addControl(child); // ✅ VBox 是容器，支持 addControl()
    _toggle->setText(_expanded ? "-" : "+");
}

void TreeNode::toggle() {
    if(_children->getNumChildren()>1){
        _expanded = !_expanded;
        _toggle->setText(_expanded ? "-" : "+");
        _children->setVisible(_expanded);
    }
}

struct App
{
    App() {
        _playing = false;
        _eci = false;
        _show._cone = false;
        _show._cov = false;
        _show._orbit = true;
        _show._label = true;
        _show._model = true;
        _show._scale = 1;
        _show._eleva = 10;
        _show._beamelev = 60;
        _speed = 10;
        readout = new ui::LabelControl();
        readout->setVertAlign(ui::Control::ALIGN_CENTER);
        setpannel = new ui::VBox();
        satlistpannel = new ui::VBox();
        _visset=true;
        _vissatlist=true;
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
    ui::CheckBoxControl* model;
    ui::VBox* setpannel;
    ui::VBox* satlistpannel;
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
    void toggelset(){_visset=!_visset;setpannel->setVisible(_visset);}
    void toggelsatlist(){_vissatlist=!_vissatlist;satlistpannel->setVisible(_vissatlist);}
    void setECI(){_eci=eci->getValue();}
    void setCone(){_show._cone=cone->getValue();}
    void setCov(){_show._cov=cov->getValue();}
    void setOrbit(){_show._orbit=orbit->getValue();}
    void setLabel(){_show._label=lab->getValue();}
    void setModel(){_show._model=model->getValue();}
    bool _eci;
    Showflags _show;
    
    bool _playing;
    int _speed;
    bool _visset;
    bool _vissatlist;
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
OE_UI_HANDLER(setModel);
OE_UI_HANDLER(toggelset);
OE_UI_HANDLER(toggelsatlist);

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
    ui::VBox* vbox = new ui::VBox();
    ui::HBox* hbox = new ui::HBox();
    hbox->addControl(new ui::ButtonControl("Set",new toggelset(app)));
    hbox->addControl(new ui::ButtonControl("SatLst",new toggelsatlist(app)));
    hbox->addControl(new ui::ButtonControl("Play", new Play(app)));
    hbox->addControl(new ui::ButtonControl("Stop", new Stop(app)));
    ui::HBox* hbox1 = new ui::HBox();
    app.lspeed=hbox1->addControl(new ui::LabelControl("speed"));
    app.speed=hbox1->addControl(new ui::HSliderControl(-100,100,10,new setSpeed(app)));
    app.speed->setWidth(150);
    hbox1->addControl(app.readout);
    vbox->addControl(hbox1);
    vbox->addControl(hbox);
    ui::Grid* grid = new ui::Grid();
    //grid->setControl(2,0,new ui::ButtonControl("Play", new Play(app)));
    //grid->setControl(3,0,new ui::ButtonControl("Stop", new Stop(app)));
    //app.lspeed=grid->setControl(4,0,new ui::LabelControl("speed"));
    //app.speed=grid->setControl(5,0,new ui::HSliderControl(-100,100,10,new setSpeed(app)));
    //app.speed->setWidth(150);
    //grid->setControl(6,0,app.readout);
    int col=0;
    app.lscale=grid->setControl(0,col,new ui::LabelControl("model scale"));
    app.scale=grid->setControl(1,col++,new ui::HSliderControl(0.1,5,1,new setScale(app)));
    app.scale->setWidth(100);
    app.lelev=grid->setControl(0,col,new ui::LabelControl("elev thrd"));
    app.elev=grid->setControl(1,col++,new ui::HSliderControl(1,70,10,new setElev(app)));
    app.elev->setWidth(100);
    app.lbeamelev=grid->setControl(0,col,new ui::LabelControl("Beam alf"));
    app.beamelev=grid->setControl(1,col++,new ui::HSliderControl(0.5,20,3,new setBeamelev(app)));
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
    grid->setControl(0,col,new ui::LabelControl("Model"));
    app.model=grid->setControl(1,col++,new ui::CheckBoxControl(true,new setModel(app)));
    
    app.setpannel->addControl(grid);
    vbox->addControl(app.setpannel);
    //vbox->addControl(app.satlistpannel);
    return vbox;
}
/*
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
*/
class ClickHandler : public osgGA::GUIEventHandler
{
public:
    //    ClickHandler(CoordinateCollector* collector,    osg::ref_ptr<PlaceNode> targetPos)
    // : _collector(collector), _targetPos(targetPos) {}
    ClickHandler(MapNode* mapnode,    osg::ref_ptr<PlaceNode> targetPos)
        : _mapNode(mapnode), _targetPos(targetPos) {}
    bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter& aa)
    {
        if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE && 
                ea.getButton() == osgGA::GUIEventAdapter::RIGHT_MOUSE_BUTTON)
        {
            osg::Vec3d world;
            if ( _mapNode->getTerrain()->getWorldCoordsUnderMouse(aa.asView(), ea.getX(), ea.getY(), world) )
            {
                GeoPoint map;
                map.fromWorld( _mapNode->getMapSRS(), world );
                std::cout << "点击位置坐标: " 
                          << "Lat=" << map.y() << ", "
                          << "Lon=" << map.x() << ", "
                          << "Alt=" << map.z() << "m" << std::endl;
                _targetPos->setPosition(map);
            }
            // 从 MouseCoordsTool 获取当前坐标
            /*const GeoPoint& point = _collector->currentPoint;
              
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
            */
            return true;
        }
        return false;
    }
    //CoordinateCollector* _collector;
    osg::ref_ptr<PlaceNode> _targetPos;
    MapNode*      _mapNode;
};
int extractLastNumber(const std::string& str) {
    int end = str.length() - 1;
    
    // 从末尾向前找到第一个非数字字符的位置
    while (end >= 0 && std::isdigit(str[end])) {
        end--;
    }
    
    // 如果没有数字，返回0或抛出异常
    if (end == str.length() - 1) {
        return 0;  // 或 throw std::invalid_argument("No digits found");
    }
    
    // 提取数字子串并转换为整数
    std::string numStr = str.substr(end + 1);
    return std::stoi(numStr);
}
std::string zero_pad(int num, int width) {
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << num;
    return oss.str();
}
void addleo(cSatellite* sat, map<string,double>& satincl, map<string,map<string,double> >& satinclraan,
            map<int,string>& satindex,int i,
            string stell,string orbmode,string satcode, double esp=0.1){
    //std::cout<<sat->Name()<<" inc "<<sat->Orbit().Inclination()*180/PI<<" raan "<<sat->Orbit().RAAN()<<" raan size "<<satinclraan.size()<<" incl size "<<satincl.size();
    double incl = sat->Orbit().Inclination();
    bool newincl = true;
    int kincl = 0;
    for(const auto &kv: satincl){
        if(kv.first.find(stell+"I")!=std::string::npos){
            kincl=extractLastNumber(kv.first);//stoi(kv.first.substr(stell.length()+1,4));
            if(fabs(incl-kv.second)<0.1){
                //std::cout<<" hit incline "<<kv.first<<" "<<kincl<<" "<<kv.second;
                newincl=false;
                break;
            }
            kincl++;
        }
    }
    if(newincl){
        //std::cout<<" new "<<kincl<<"  incline "<<incl<<" I"<<to_string(kincl);
        satincl[stell+"I"+to_string(kincl)]=incl;
    }
    double raan = sat->Orbit().RAAN();
    bool newraan=true;
    int k=0;
    int ki=0;
    for(const auto &kv: satinclraan[stell+"I"+to_string(kincl)]){
        k=extractLastNumber(kv.first);//stoi(kv.first.substr(stell.length()+1,4));
        if(fabs(raan-kv.second)<esp){
            //std::cout<<" hit "<<kv.first<<" "<<k<<" "<<kv.second<<std::endl;
            newraan=false;
            break;
        }
        k++;
    }
    
    if(newraan){
        //std::cout<<" new "<<k<<"  raan "<<raan<<" P"<<to_string(k)<<std::endl;
        satinclraan[stell+"I"+to_string(kincl)]["P"+zero_pad(k,2)]=raan;
    }
    satindex[i]=stell+"I"+to_string(kincl)+"P"+zero_pad(k,2);
}
void addmeo(cSatellite* sat,map<string,TreeNode*>& stellnode, map<string,double>& satraan,SatelliteObj* satobj,int i,
            string stell,string orbmode,string satcode,string stype,double esp=0.5){
    std::cout<<sat->Name()<<" inc "<<sat->Orbit().Inclination()<<" raan "<<sat->Orbit().RAAN()<<" size "<<satraan.size();
    double raan = sat->Orbit().RAAN();
    bool newraan=true;
    int k=0;
    int kmax=k;
    for(const auto & kv: satraan){
        //std::cout<< kv.first<<" -> "<<kv.second<<"  "<<kv.first.substr(stell.length()+1,4)<<std::endl;
        if(kv.first.find(stell+"P")!=std::string::npos){
            k=stoi(kv.first.substr(stell.length()+1,4));
            if(fabs(raan-kv.second)<esp){
                std::cout<<" hit "<<kv.first<<" "<<k<<" "<<kv.second<<std::endl;
                newraan=false;
                break;
            }
            k++;
            kmax=kmax>k?kmax:k;
        }
    }
    if(newraan){
        std::cout<<" new "<<kmax<<"  raan "<<raan<<" P"<<to_string(kmax)<<std::endl;
        k=kmax;
        satraan[stell+"P"+to_string(k)]=raan;
        stellnode[stell+"P"+to_string(k)]=new TreeNode(stype+"P"+to_string(k));
        stellnode[stell]->addChild(stellnode[stell+"P"+to_string(k)]);
    }
    stellnode[stell+"P"+to_string(k)]->addChild(new TreeNode(orbmode+satcode,i,satobj));
    
}
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
        
        //osgEarth::SimpleSky::SimpleSkyOptions skyOptions;
        //osgEarth::SimpleSky::SimpleSkyNode* skyNode = new osgEarth::SimpleSky::SimpleSkyNode(skyOptions);
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
        //container->addChild(createUI(app));
        ui::VBox* vbox = new VBox();
        vbox->addControl(createUI(app));
        VBox* satroot = new VBox();
        //TreeNode* nodebd3 = new TreeNode("BEIDOU-3");
        //satroot->addControl(nodebd3);
        app.satlistpannel->addControl(satroot);
        vbox->addControl(app.satlistpannel);
        container->addChild(vbox);
        
        viewer.addEventHandler(new ClickHandler(mapNode,app.targetPos));
        // A lat/long SRS for specifying points.
        //const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
        //std::map<int,osg::ref_ptr<osg::Geode> > satobj;
        osg::ref_ptr<osg::Geode> geode = new osg::Geode();
        std::map<int,SatelliteObj* > satobjmap;
        const DateTime& t = app.sky->getDateTime();
        CelestialBody sun = ephemeris->getSunPosition(t);
        osg::Vec3d sunPos = sun.geocentric;
        sunPos.normalize();
        //TreeNode* nodebd3geo = new TreeNode("GEO");
        //nodebd3->addChild(nodebd3geo);
        //TreeNode* nodebd3igso = new TreeNode("IGSO");
        //nodebd3->addChild(nodebd3igso);
        //TreeNode* nodebd3meo = new TreeNode("MEO");
        //nodebd3->addChild(nodebd3meo);
        map<string,TreeNode*> stellnode;
        map<string,double> satraan;
        map<string,double> satincl;
        map<string,map<string,double> > satinclraan;
        map<int,string> satindex;
        //vector<TreeNode*> satplane;
        int planenum=0;
        for(int i=0;i<satlist.size();i++)
        {
            cSatellite* sat = satlist[i];
            string stell,orbmode,satcode;
            std::regex pattern("([^\\s]+)\\s+([^\\s]+)\\s+(\\([^\\a]+\\))");//BEIDOU COSMOS GPS
            std::smatch match;
            std::string input = sat->Name();
            //std::cout<<input<<":";
            if(std::regex_search(input,match,pattern)){
                //std::cout<<"regex:"<<input<<" res "<<match[1]<<" "<<match[2]<<" "<<match[3]<<std::endl;
                stell=match[1];
                orbmode=match[2];
                satcode=match[3];
            }else{
                pattern="([^\\s]+)\\s+(\\((\\w+)[^\\a]+\\))";//GALLILEO QZSS NVS
                if(std::regex_search(input,match,pattern)){
                    stell = match[3];
                    orbmode = match[1];
                    satcode = match[2];
                    if(stell=="NVS") stell="IRNSS";
                    else if(stell=="DARKSA") stell="STARLINK";
                }else{
                    pattern="((\\w+)-\\w+)";
                    if(std::regex_search(input,match,pattern)){
                        stell = match[2];
                        orbmode = match[1];
                        satcode = "";
                    }
                }
            }
            //std::cout<<"regex: "<<stell<<"   |   "<<orbmode<<"   |   "<<satcode<<std::endl;
            if(stell=="BEIDOU-2"||stell=="BEIDOU-3")
            {
                if(stellnode[stell]==nullptr){
                    stellnode[stell]=new TreeNode(stell);
                    satroot->addControl(stellnode[stell]);
                    stellnode[stell+"GEO"]=new TreeNode("GEO");
                    stellnode[stell]->addChild(stellnode[stell+"GEO"]);
                    stellnode[stell+"IGSO"]=new TreeNode("IGSO");
                    stellnode[stell]->addChild(stellnode[stell+"IGSO"]);
                    stellnode[stell+"MEO"]=new TreeNode("MEO");
                    stellnode[stell]->addChild(stellnode[stell+"MEO"]);
                }
                SatelliteObj * satobj = new SatelliteObj(sat,mapNode,geode,t,sunPos,app.targetPos);
                root->addChild(satobj->pat);
                satobjmap[i]=satobj;
                if(orbmode.find("I")!=std::string::npos)
                {
                    if(stell=="BEIDOU-2")
                        addmeo(sat,stellnode,satraan,satobj,i,stell+"IGSO",orbmode,satcode,"IGSO");
                    else
                        stellnode[stell+"IGSO"]->addChild(new TreeNode(orbmode+satcode,i,satobj));
                }else if(orbmode.find("G")!=std::string::npos){
                    stellnode[stell+"GEO"]->addChild(new TreeNode(orbmode+satcode,i,satobj));
                }else{
                    addmeo(sat,stellnode,satraan,satobj,i,stell+"MEO",orbmode,satcode,"MEO");
                }
                
            }else if(stell=="COSMOS"||stell=="GPS"||stell=="GALILEO"){
                if(stellnode[stell+"MEO"]==nullptr){
                    stellnode[stell+"MEO"]=new TreeNode(stell);
                    satroot->addControl(stellnode[stell+"MEO"]);
                }
                SatelliteObj * satobj = new SatelliteObj(sat,mapNode,geode,t,sunPos,app.targetPos);
                root->addChild(satobj->pat);
                addmeo(sat,stellnode,satraan,satobj,i,stell+"MEO",orbmode,satcode,"MEO");
                satobjmap[i]=satobj;
            }else if(stell=="QZSS"||stell=="IRNSS"){
                if(stellnode[stell+"GSO"]==nullptr){
                    stellnode[stell+"GSO"]=new TreeNode(stell);
                    satroot->addControl(stellnode[stell+"GSO"]);
                }
                SatelliteObj * satobj = new SatelliteObj(sat,mapNode,geode,t,sunPos,app.targetPos);
                root->addChild(satobj->pat);
                if(stell=="QZSS")
                    stellnode[stell+"GSO"]->addChild(new TreeNode(orbmode+satcode,i,satobj));
                else
                    addmeo(sat,stellnode,satraan,satobj,i,stell+"GSO",orbmode,satcode,"GSO");
                satobjmap[i]=satobj;
                //}else if(stell=="STARLINK"){
                
            }else if(stell=="STARLINK"){
                addleo(sat,satincl,satinclraan,satindex,i,stell,orbmode,satcode);
            }
            else{
                std::cout<<stell<<" "<<orbmode<<" "<<satcode<<std::endl;
                stell="Other";
                if(satcode.empty()) satcode=sat->Name();
                if(stellnode[stell+"GSO"]==nullptr){
                    stellnode[stell+"GSO"]=new TreeNode(stell);
                    satroot->addControl(stellnode[stell+"GSO"]);
                }
                SatelliteObj * satobj = new SatelliteObj(sat,mapNode,geode,t,sunPos,app.targetPos);
                root->addChild(satobj->pat);
                stellnode[stell+"GSO"]->addChild(new TreeNode(orbmode+satcode,i,satobj));
                satobjmap[i]=satobj;
                
            }
        }
        if(satindex.size()>0){
            map<string,string> igp;
            string stell="STARLINK";
            stellnode[stell]=new TreeNode(stell);
            satroot->addControl(stellnode[stell]);
            cout<<"stellnode "<<stellnode.size()<<endl;
            for(const auto & kv: satinclraan){
                //cout<<kv.first<<" "<<kv.second.size()<<endl;
                stellnode[kv.first]=new TreeNode(kv.first.substr(8,2)+"/"+to_string((int)(satincl[kv.first]*180/PI)));
                stellnode[stell]->addChild(stellnode[kv.first]);
                std::map<double,string> mraan;
                for(const auto & kv1:kv.second){
                    mraan[kv1.second]=kv1.first;
                }
                int g=0;
                int p=0;
                if(kv.second.size()){
                    int c=0;
                    for(const auto & kv1:mraan){
                        if(c%10==0){
                            g++;
                            stellnode[kv.first+"G"+to_string(g)]=new TreeNode("G"+to_string(g));
                            stellnode[kv.first]->addChild(stellnode[kv.first+"G"+to_string(g)]);
                            p=0;};
                        //cout<<kv.first<<" "<<kv1.second<<" g"<<g<<" p"<<c<<" "<<kv1.first*180/PI<<endl;
                        igp[kv.first+kv1.second]=kv.first+"G"+to_string(g)+"P"+zero_pad(c,2);
                        stellnode[kv.first+"G"+to_string(g)+"P"+zero_pad(c,2)]=new TreeNode("P"+to_string(c));
                        stellnode[kv.first+"G"+to_string(g)]->addChild(stellnode[kv.first+"G"+to_string(g)+"P"+zero_pad(c,2)]);
                        c++;
                    }
                    
                }
            }
            
            //std::cout<<"satindex size "<<satindex.size()<<std::endl;
            for(const auto & kv:satindex){
                //	int i=stoi(igp[kv.second].substr(9,1));
                //	int g=stoi(igp[kv.second].substr(11,1));
                //	int p=stoi(igp[kv.second].substr(13,2));
                //        cout<<kv.first<<" "<<kv.second<<" "<<igp[kv.second]<<" "<<i<<" "<<g<<" "<<p<<endl;
                SatelliteObj * satobj = new SatelliteObj(satlist[kv.first],mapNode,geode,t,sunPos,app.targetPos);
                root->addChild(satobj->pat);
                stellnode[igp[kv.second]]->addChild(new TreeNode(satlist[kv.first]->Name(),kv.first,satobj));
                satobjmap[kv.first]=satobj;
            }
        }
        root->addChild(geode);
        
        /*
for (const auto& pair : satobjmap)
{
  std::cout << "Key: " << pair.first << ", Value: " << pair.second->sat->Name() << std::endl;
}*/
        
        
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
