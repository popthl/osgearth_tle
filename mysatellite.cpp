#include "mysatellite.h"
#include <osgEarth/LabelNode>
SatelliteObj::SatelliteObj(cSatellite* sattle,MapNode* mapnode,osg::ref_ptr<osg::Geode> geode,time_t t0){
    mapNode=mapnode;
    sat=sattle;
    satgroup=geode;
    boxTransform = new osg::MatrixTransform;

    const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
    double Tcycle = 2*PI/sat->Orbit().MeanMotion();
    double dt = Tcycle/3*60;
    //printf("Tcycle %f dt %f t0 %ld\n",Tcycle*60,dt,t0);
    vector<osg::Vec3d> pt;
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    orbit = new osg::Geometry;
    osg::Vec3d llh0,xyz;
    for(int i=0;i<3;i++)
    {
        //printf("t %d %s",t0,ctime(&t0));
        vector<double> llh = SatEcf(*sat,t0);
        double dl = i*dt*OMEGA_E;//IAU_EARTH_ANGULAR_VELOCITY;
        //cout<<llh[0]<<" "<<llh[1]<<" "<<llh[2]<<" "<<dl<<endl;
	llh0 =osg::Vec3d(llh[0]+dl, llh[1], llh[2]*1000);
	pt.push_back(llh0);
	if(i>0){
	    Vec3dVector out;
	    TessellateOperator::tessellateGeo(pt[i-1],pt[i],20,GEOINTERP_GREAT_CIRCLE,out);
	    for(int k=0;k<out.size();k++){
		geoSRS->transformToWorld(out[k], xyz);
		vertices->push_back(xyz);
		//printf("%d : %lf %lf %lf\n",k,out[k].x(),out[k].y(),out[k].z());
	    }
	}
        t0+=dt;
    }
    Vec3dVector out;
    //printf("pt0: %lf %lf %lf\n",pt[0].x(),pt[0].y(),pt[0].z());
    TessellateOperator::tessellateGeo(pt[2],pt[0],20,GEOINTERP_GREAT_CIRCLE,out);
    for(int i=0;i<out.size();i++){
	geoSRS->transformToWorld(out[i], xyz);
	vertices->push_back(xyz);
	//printf("%d : %lf %lf %lf\n",i,out[i].x(),out[i].y(),out[i].z());
    }
    orbit->setVertexArray(vertices);
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    orbit->setColorArray(colors);
    orbit->setColorBinding(osg::Geometry::BIND_OVERALL);

    // 指定图元类型（这里是 LINE_LOOP）
    orbit->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, vertices->size()));

    // 设置渲染状态
    osg::ref_ptr<osg::StateSet> stateSet = orbit->getOrCreateStateSet();
    osg::ref_ptr<osg::LineWidth> lineWidth = new osg::LineWidth;
    lineWidth->setWidth(2.0f);
    stateSet->setAttributeAndModes(lineWidth, osg::StateAttribute::ON);

    boxwing = createColoredBox(osg::Vec3(0,0,0),osg::Vec3(RADIUS,RADIUS*1,RADIUS*2));

    Style labelStyle;
    labelStyle.getOrCreate<TextSymbol>()->alignment() = TextSymbol::ALIGN_CENTER_CENTER;
    labelStyle.getOrCreate<TextSymbol>()->fill()->color() = Color::Yellow;
    label = new LabelNode(sat->Name(), labelStyle);
    //label->setPosition();
    label->setPosition(GeoPoint(geoSRS, pt[0].x(), pt[0].y(),pt[0].z()));

    geode->addDrawable(orbit);
    //geode->addDrawable(boxwing);
    mapNode->addChild(label);
    osg::Matrix translationMatrix = osg::Matrix::translate((*vertices)[0]);
    boxTransform->setMatrix(translationMatrix);
    boxTransform->addChild(boxwing);

}
void SatelliteObj::setposition(time_t t0){
    //osg::Vec3 lastpos =(* dynamic_cast<osg::Vec3Array*>(orbit->getVertexArray()) )[0];
    const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
    double Tcycle = 2*PI/sat->Orbit().MeanMotion();
    double dt = Tcycle/3*60;
    //printf("Tcycle %f dt %f t0 %ld\n",Tcycle*60,dt,t0);
    vector<osg::Vec3d> pt;
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    osg::Vec3d llh0,xyz;
    for(int i=0;i<3;i++)
    {
        //printf("t %d %s",t0,ctime(&t0));
        vector<double> llh = SatEcf(*sat,t0);
        double dl = i*dt*OMEGA_E;//IAU_EARTH_ANGULAR_VELOCITY;
        //cout<<llh[0]<<" "<<llh[1]<<" "<<llh[2]<<" "<<dl<<endl;
	llh0 =osg::Vec3d(llh[0]+dl, llh[1], llh[2]*1000);
	pt.push_back(llh0);
	if(i>0){
	    Vec3dVector out;
	    TessellateOperator::tessellateGeo(pt[i-1],pt[i],20,GEOINTERP_GREAT_CIRCLE,out);
	    for(int k=0;k<out.size();k++){
		geoSRS->transformToWorld(out[k], xyz);
		vertices->push_back(xyz);
		//printf("%d : %lf %lf %lf\n",k,out[k].x(),out[k].y(),out[k].z());
	    }
	}
        t0+=dt;
    }
    Vec3dVector out;
    //printf("pt0: %lf %lf %lf\n",pt[0].x(),pt[0].y(),pt[0].z());
    TessellateOperator::tessellateGeo(pt[2],pt[0],20,GEOINTERP_GREAT_CIRCLE,out);
    for(int i=0;i<out.size();i++){
	geoSRS->transformToWorld(out[i], xyz);
	vertices->push_back(xyz);
	//printf("%d : %lf %lf %lf\n",i,out[i].x(),out[i].y(),out[i].z());
    }
    orbit->setVertexArray(vertices);
    label->setPosition(GeoPoint(geoSRS, pt[0].x(), pt[0].y(),pt[0].z()));
    osg::Matrix translationMatrix = osg::Matrix::translate((*vertices)[0]);
    boxTransform->setMatrix(translationMatrix);

}
/*osg::ref_ptr<osg::Geode> createsatellite(cSatellite* sat,MapNode* mapNode,time_t t0){
    const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
    osg::ref_ptr<osg::Geode> geode = new osg::Geode();
    double Tcycle = 2*PI/sat->Orbit().MeanMotion();
    double dt = Tcycle/3*60;
    //printf("Tcycle %f dt %f t0 %ld\n",Tcycle*60,dt,t0);
    vector<osg::Vec3d> pt;
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    osg::ref_ptr<osg::Geometry> orbit = new osg::Geometry;
    osg::Vec3d llh0,xyz;
    for(int i=0;i<3;i++)
    {
        //printf("t %d %s",t0,ctime(&t0));
        vector<double> llh = SatEcf(*sat,t0);
        double dl = i*dt*OMEGA_E;//IAU_EARTH_ANGULAR_VELOCITY;
        //cout<<llh[0]<<" "<<llh[1]<<" "<<llh[2]<<" "<<dl<<endl;
	llh0 =osg::Vec3d(llh[0]+dl, llh[1], llh[2]*1000);
	pt.push_back(llh0);
	if(i>0){
	    Vec3dVector out;
	    TessellateOperator::tessellateGeo(pt[i-1],pt[i],20,GEOINTERP_GREAT_CIRCLE,out);
	    for(int k=0;k<out.size();k++){
		geoSRS->transformToWorld(out[k], xyz);
		vertices->push_back(xyz);
		//printf("%d : %lf %lf %lf\n",k,out[k].x(),out[k].y(),out[k].z());
	    }
	}
        t0+=dt;
    }
    Vec3dVector out;
    //printf("pt0: %lf %lf %lf\n",pt[0].x(),pt[0].y(),pt[0].z());
    TessellateOperator::tessellateGeo(pt[2],pt[0],20,GEOINTERP_GREAT_CIRCLE,out);
    for(int i=0;i<out.size();i++){
	geoSRS->transformToWorld(out[i], xyz);
	vertices->push_back(xyz);
	//printf("%d : %lf %lf %lf\n",i,out[i].x(),out[i].y(),out[i].z());
    }
    osg::ref_ptr<osg::Geometry> coloredBox = createColoredBox((*vertices)[0],osg::Vec3(RADIUS,RADIUS*1,RADIUS*2));
    geode->addDrawable(coloredBox);
    //mapNode->addChild(coloredBox);
    // Style our labels:
    Style labelStyle;
    labelStyle.getOrCreate<TextSymbol>()->alignment() = TextSymbol::ALIGN_CENTER_CENTER;
    labelStyle.getOrCreate<TextSymbol>()->fill()->color() = Color::Yellow;
    LabelNode* label = new LabelNode(sat->Name(), labelStyle);
    //label->setPosition();
    label->setPosition(GeoPoint(geoSRS, pt[0].x(), pt[0].y(),pt[0].z()));
    mapNode->addChild(label);


    /*osg::ref_ptr<osgText::Text> text = new osgText::Text;
    text->setCharacterSize(20000.0);
    text->setPosition(xyz);
    text->setText(sat->Name());
    mapNode->addChild(text);*

    orbit->setVertexArray(vertices);
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    orbit->setColorArray(colors);
    orbit->setColorBinding(osg::Geometry::BIND_OVERALL);

    // 指定图元类型（这里是 LINE_LOOP）
    orbit->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, vertices->size()));

    // 设置渲染状态
    osg::ref_ptr<osg::StateSet> stateSet = orbit->getOrCreateStateSet();
    osg::ref_ptr<osg::LineWidth> lineWidth = new osg::LineWidth;
    lineWidth->setWidth(2.0f);
    stateSet->setAttributeAndModes(lineWidth, osg::StateAttribute::ON);
    geode->addDrawable(orbit);

    return geode;

}*/
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
                }else{*/
                    linename=line1.substr(1,7);
                }
                //std::cout<<noradid<<std::endl;
                cTle tle(linename, line1, line2);
                cSatellite* sat=new cSatellite(tle);
                time_t ep = sat->Orbit().Epoch().ToTime();
                //printf("sat %s ep %s",sat->Name().c_str(),ctime(&ep));
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
