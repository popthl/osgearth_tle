#include "mysatellite.h"
#include <osgEarth/LabelNode>
#include <osg/BlendFunc>
#include <osgEarth/FeatureNode>
#include <osgEarth/GeoData>

osg::ref_ptr<osg::MatrixTransform> createRotatingSail(double angle);
osg::ref_ptr<osg::Node> createColoredCube();
vector<vector<double> > SatBeam(vector<double>spos, double alf, double tlon, double tlat,double eleva);

SatelliteObj::SatelliteObj(cSatellite* sattle,MapNode* mapnode,osg::ref_ptr<osg::Geode> geode,DateTime t0,osg::Vec3d sunPos,osg::ref_ptr<PlaceNode> targetPos){
    mapNode=mapnode;
    sat=sattle;
    satgroup=geode;
    target = targetPos;
    _vis=true;
    //boxTransform = new osg::MatrixTransform;
    pat = new osg::PositionAttitudeTransform();
    showflag._scale = 0.5;
    showflag._eleva = 35;//deg
    showflag._beamelev = 75;//deg
    box = createColoredCube();
    wing = createRotatingSail(1.0);

    Style labelStyle;
    labelStyle.getOrCreate<TextSymbol>()->alignment() = TextSymbol::ALIGN_CENTER_CENTER;
    labelStyle.getOrCreate<TextSymbol>()->fill()->color() = Color::Yellow;
    label = new LabelNode(sat->Name(), labelStyle);

    orbit = new osg::Geometry;
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));
    orbit->setColorArray(colors);
    orbit->setColorBinding(osg::Geometry::BIND_OVERALL);
    osg::ref_ptr<osg::StateSet> stateSet = orbit->getOrCreateStateSet();
    osg::ref_ptr<osg::LineWidth> lineWidth = new osg::LineWidth;
    lineWidth->setWidth(2.0f);
    stateSet->setAttributeAndModes(lineWidth, osg::StateAttribute::ON);

    cov = new osg::Geometry;
    osg::ref_ptr<osg::Vec4Array> colors1 = new osg::Vec4Array;
    colors1->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 0.5f));
    cov->setColorArray(colors1);
    cov->setColorBinding(osg::Geometry::BIND_OVERALL);
    osg::ref_ptr<osg::StateSet> stateSet1 = cov->getOrCreateStateSet();
    //stateSet1->setAttributeAndModes(lineWidth, osg::StateAttribute::ON);

    stateSet1->setMode(GL_BLEND, osg::StateAttribute::ON);
    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc;
    blendFunc->setFunction(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    stateSet1->setAttributeAndModes(blendFunc.get(), osg::StateAttribute::ON);

    osgEarth::Geometry* geom = new osgEarth::Polygon();
    osgEarth::Style style;
    style.getOrCreate<PolygonSymbol>()->fill()->color() = Color(Color::White, 0.5);
    osgEarth::Feature* feature = new osgEarth::Feature(geom, mapNode->getMapSRS()->getGeographicSRS());
    featureNode = new osgEarth::FeatureNode(feature,style);
    mapNode->addChild(featureNode);
    //geode->addDrawable(orbit);
    showflag._orbit=true;
    showflag._cone = true;
    showflag._cov = true;
    showflag._label = true;
    setposition(t0,sunPos);

//printf("ori satgroup %d\n",satgroup->getNumDrawables());
    mapNode->addChild(label);
    pat->addChild(box);
    pat->addChild(wing);
}
void SatelliteObj::setvisible(bool vis){
    _vis=vis;
    if(_vis){
	cov->setNodeMask(0xffffffff);
	orbit->setNodeMask(0xffffffff);
	box->setNodeMask(0xffffffff);
	label->setNodeMask(0xffffffff);
	wing->setNodeMask(0xffffffff);
	featureNode->setNodeMask(0xffffffff);
    }
    else{
	cov->setNodeMask(0);
	orbit->setNodeMask(0);
	box->setNodeMask(0);
	label->setNodeMask(0);
	wing->setNodeMask(0);
	featureNode->setNodeMask(0);
    }
}
int SatelliteObj::setposition(DateTime t0,osg::Vec3d sunPos){
if(!_vis)return 0;
    const SpatialReference* geoSRS = mapNode->getMapSRS()->getGeographicSRS();
    double Tcycle = 2*PI/sat->Orbit().MeanMotion()/60;
osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
osg::Vec3d llh0,xyz;
DateTime t=t0;
double h0;
cJulian cjt0(t0.asTimeStamp());
	vector<osg::Vec3d> pt;
{
	    cEciTime sateci = sat->PositionEci(cjt0);
	    cEci pos = sateci;
	    cGeo satgeo(pos, cjt0);
	    llh0=osg::Vec3d(satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm()*1000);
            geoSRS->transformToWorld(llh0, xyz);
            vertices->push_back(xyz);
	    h0=satgeo.AltitudeKm();
	    pt.push_back(llh0);
}
//printf("showflag %d %d %d, vertices %d satgroup %d\n",showflag._orbit,showflag._cone,showflag._cov,
//    vertices->size(),satgroup->getNumDrawables());
if(showflag._orbit){
    //printf("Tcycle %f\n",Tcycle*3600);
    //cout<<"t0: "<<t0.asCompactISO8601()<<endl;
    int npoint = 3;
    double ra = sat->Orbit().Apogee();
    double rp = sat->Orbit().Perigee();

    if(sat->Orbit().Eccentricity()>0.02){
	npoint=40;
	while(t<t0+Tcycle){
	    double s = 1+5*(ra-h0)/(ra-rp);
	    //printf("s %f\n",s);
	    t=t+Tcycle/npoint/s;
	    cJulian cjt(t.asTimeStamp());
	    cEciTime sateci = sat->PositionEci(cjt);
	    cEci pos = sateci;
	    cGeo satgeo(pos, cjt0);
	    //printf("llh: %lf %lf %lf\n",satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm());
	    llh0=osg::Vec3d(satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm()*1000);
            geoSRS->transformToWorld(llh0, xyz);
            vertices->push_back(xyz);
	    h0=satgeo.AltitudeKm();
	}
    }
    else{
	npoint=3;

	double dt = Tcycle/npoint;
        //printf("Tcycle %f dt %f \n",Tcycle*60,dt);

    	for(int i=1;i<npoint;i++)
    	{
	    t=t0+i*dt;
	    //cout<<t.asCompactISO8601()<<endl;
	    cJulian cjt(t.asTimeStamp());
	    cEciTime sateci = sat->PositionEci(cjt);
	    cEci pos = sateci;
	    cGeo satgeo(pos, cjt0);
            //printf("llh: %lf %lf %lf\n",satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm());
	    llh0=osg::Vec3d(satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm()*1000);
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
	}
	Vec3dVector out;
	TessellateOperator::tessellateGeo(pt[npoint-1],pt[0],20,GEOINTERP_GREAT_CIRCLE,out);
	for(int i=0;i<out.size();i++){
	    geoSRS->transformToWorld(out[i], xyz);
	    vertices->push_back(xyz);
	    //printf("%d : %lf %lf %lf\n",i,out[i].x(),out[i].y(),out[i].z());
	}
    }
    //printf("vertices %d\n",vertices->size());
    if(orbit->getNumPrimitiveSets()==0)
	orbit->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP, 0, vertices->size()));
    orbit->setVertexArray(vertices);
    if(!satgroup->containsDrawable(orbit))
	satgroup->addDrawable(orbit);
}else{
	satgroup->removeDrawable(orbit);
	    t=t0+Tcycle/6;
	    cJulian cjt(t.asTimeStamp());
	    cEciTime sateci = sat->PositionEci(cjt);
	    cEci pos = sateci;
	    cGeo satgeo(pos, cjt0);
	    llh0=osg::Vec3d(satgeo.LongitudeDeg(),satgeo.LatitudeDeg(),satgeo.AltitudeKm()*1000);
            geoSRS->transformToWorld(llh0, xyz);
            vertices->push_back(xyz);

}
    GeoPoint pos0;
    pos0.fromWorld(geoSRS,(*vertices)[0]);
    label->setPosition(pos0);
if(showflag._label)
	label->setNodeMask(0xffffffff);//Visible(showflag._label);
else
	label->setNodeMask(0);
    setatt(vertices,sunPos);
//printf("showflag %d %d %d, vertices %d satgroup %d\n",showflag._orbit,showflag._cone,showflag._cov,
//    vertices->size(),satgroup->getNumDrawables());

//-------------------add satbeam
if(showflag._cone || showflag._cov){
    vector<double>spos;
    spos.push_back(pos0.x());
    spos.push_back(pos0.y());
    spos.push_back(pos0.z());
    double alf = asin(6378140/pos0.z()*cos(showflag._beamelev*PI/180));
//printf("llh: %lf %lf %lf\n",spos[0],spos[1],spos[2]);
    bool icov;
    GeoPoint tpos = target->getPosition();
    vector<vector<double> > covpt = SatBeam(spos, alf, tpos.x(),tpos.y(),showflag._eleva);//pos0.x(), pos0.y());
//printf("covpt %d\n",covpt.size());
    osgEarth::Geometry* geom = new osgEarth::Polygon();
    if(covpt.size()>0){
	osg::ref_ptr<osg::Vec3Array> vertices1 = new osg::Vec3Array;
	vertices1->push_back((*vertices)[0]);
	for(int i=0;i<covpt.size();i++)
	{
	    osg::Vec3d pt(covpt[i][0],covpt[i][1],100000);
	    geoSRS->transformToWorld(pt, xyz);
	    vertices1->push_back(xyz);
//printf("%d : %lf %lf %lf  %lf %lf %lf\n",i,xyz.x(),xyz.y(),xyz.z(),pt.x(),pt.y(),pt.z());
	    if(showflag._cov)
		geom->push_back(pt);
	}
	vertices1->push_back((*vertices1)[1]);
    	if(cov->getNumPrimitiveSets()==0)
	    cov->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::TRIANGLE_FAN, 0, vertices1->size()));
    	cov->setVertexArray(vertices1);
	icov = true;
    }else{
	icov=false;
    }
    if(showflag._cone){
	if(icov){
	    if(!satgroup->containsDrawable(cov))
		satgroup->addDrawable(cov);
	}else
		satgroup->removeDrawable(cov);
    }else
		satgroup->removeDrawable(cov);
    osgEarth::Feature* feature = new osgEarth::Feature(geom,geoSRS);
    featureNode->setFeature(feature);
}else{
 if(showflag._cone==false){
	satgroup->removeDrawable(cov);
 }
 if(showflag._cov==false){
    osgEarth::Geometry* geom = new osgEarth::Polygon();
    osgEarth::Feature* feature = new osgEarth::Feature(geom,geoSRS);
    featureNode->setFeature(feature);
 }
}
//------------------------
    return vertices->size();
}
void SatelliteObj::setatt(osg::ref_ptr<osg::Vec3Array> vertices,osg::Vec3d sunPos)
{
    /*printf("sunPos: %lf,%lf,%lf\n",sunPos.x(),sunPos.y(),sunPos.z());
    printf("satPos: %lf,%lf,%lf\n",(*vertices)[0].x()/(*vertices)[0].length(),(*vertices)[0].y()/(*vertices)[0].length(),(*vertices)[0].z()/(*vertices)[0].length());
    printf("satPos: %lf,%lf,%lf\n",(*vertices)[10].x()/(*vertices)[10].length(),(*vertices)[10].y()/(*vertices)[10].length(),(*vertices)[10].z()/(*vertices)[10].length());
*/
   // rotate satellite body axis z, point to earth center
   double dot=(*vertices)[0] * osg::Vec3(0,0,1);
   double r=(*vertices)[0].length();
   double angle = PI-std::acos(dot/r);
   pat->setPosition((*vertices)[0]);
   osg::Vec3 axisz=(*vertices)[0]^osg::Vec3(0,0,1);
   osg::Quat rotation(angle, axisz);
   //pat->setAttitude(rotation);
   osg::Vec3 axisy0 = rotation*osg::Vec3(0,1,0);
   osg::Vec3 axisx0 = rotation*osg::Vec3(1,0,0);
   int id=1;
   if(showflag._orbit)
	id=10;
   osg::Vec3 axisy=(*vertices)[0]^(*vertices)[id];//orbit normal direction
   axisy.normalize();
//    printf("axisy: %lf,%lf,%lf\n",axisy.x(),axisy.y(),axisy.z());
   dot=axisy*axisy0;
   double dot1 = axisy*axisx0;
   angle = std::acos(dot);
   if(std::isnan(angle))
	angle=PI;
   //printf("angle: %f,dot %f dot1:%f\n",angle*180/PI,dot,dot1);
   if((dot>0 && dot1<0)||(dot<0 && dot1<0))
	angle=-angle;
   //body x align to move dirction

   //compute yaw and sail pan rotate angle in yaw steering mode
   osg::Vec3 orbity=axisy^sunPos;
   orbity.normalize();
//    printf("orbity: %lf,%lf,%lf\n",orbity.x(),orbity.y(),orbity.z());
   double sinu = (*vertices)[0] * orbity/r;
   double beta = asin(axisy*sunPos);
   double yaw = atan2(tan(beta),sinu);
   angle-=yaw;
   osg::Quat rotatey(angle,(*vertices)[0]);
   pat->setAttitude(rotation*rotatey);
   pat->setScale(osg::Vec3(showflag._scale*RADIUS, showflag._scale*RADIUS, showflag._scale*RADIUS));
   double anglewing=acos(cos(beta)*sqrt(1-sinu*sinu));
    dot = sunPos*(*vertices)[0];
   if(dot>0)
       anglewing=-anglewing;
//printf("beta %f, yaw %f, sinu %f, anglewing %f\n",beta*180/PI,yaw*180/PI,sinu,anglewing*180/PI);
   wing->setMatrix(osg::Matrix::rotate(anglewing, 0, 1, 0));
}

#define pi 3.14159265358979323846
#define undefined 999999.1
/* -----------------------------------------------------------------------------
*
*                           function dot
*
*  this function finds the dot product of two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    dot         - result
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

    double  dot
        (
        double x[3], double y[3]
        )
    {
        return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
    }  // dot

/* -----------------------------------------------------------------------------
*
*                           function mag
*
*  this procedure finds the magnitude of a vector.  the tolerance is set to
*    0.000001, thus the 1.0e-12 for the squared test of underflows.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec       - vector
*
*  outputs       :
*    vec       - answer stored in function return
*
*  locals        :
*    none.
*
*  coupling      :
*    none.
* --------------------------------------------------------------------------- */

    double  mag
        (
        double x[3]
        )
    {
        return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    }  // mag

/* -----------------------------------------------------------------------------
*
*                           procedure cross
*
*  this procedure crosses two vectors.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec1        - vector number 1
*    vec2        - vector number 2
*
*  outputs       :
*    outvec      - vector result of a x b
*
*  locals        :
*    none.
*
*  coupling      :
*    none
 ---------------------------------------------------------------------------- */

    void    cross
        (
        double vec1[3], double vec2[3], double outvec[3]
        )
    {
        outvec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
        outvec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
        outvec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    }  // cross

/* -----------------------------------------------------------------------------
*
*                           procedure norm
*
*  this procedure calculates a unit vector given the original vector.  if a
*    zero vector is input, the vector is set to zero.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs          description                    range / units
*    vec         - vector
*
*  outputs       :
*    outvec      - unit vector
*
*  locals        :
*    i           - index
*    small       - value defining a small value
*    magv        - magnitude of the vector
*
*  coupling      :
*    mag           magnitude of a vector
* --------------------------------------------------------------------------- */

    void    norm
        (
        double vec[3],
        double outvec[3]
        )
    {
        const double small = 0.000001;
        double magv;
        int i;

        magv = mag(vec);
        if (magv > small)
        {
            for (i = 0; i <= 2; i++)
                outvec[i] = vec[i] / magv;
        }
        else
        for (i = 0; i <= 2; i++)
            outvec[i] = 0.0;
    }  // norm
    /* -----------------------------------------------------------------------------
    *
    *                           procedure addvec
    *
    *  this procedure adds two vectors possibly multiplied by a constant.
    *
    *  author        : david vallado                  719-573-2600    1 mar 2001
    *
    *  inputs          description                    range / units
    *    a1          - constant multiplier
    *    a2          - constant multiplier
    *    vec1        - vector number 1
    *    vec2        - vector number 2
    *
    *  outputs       :
    *    outvec      - vector result of a + b
    *
    *  locals        :
    *    row         - index
    *
    *  coupling      :
    *     none
    * --------------------------------------------------------------------------- */

        void    addvec
            (
            double a1, double vec1[3],
            double a2, double vec2[3],
            double vec3[3]
            )
        {
            int row;

            for (row = 0; row <= 2; row++)
            {
                vec3[row] = 0.0;
                vec3[row] = a1* vec1[row] + a2* vec2[row];
            }
        }  // addvec

        /* -----------------------------------------------------------------------------
        *
        *                           procedure matvecmult
        *
        *  this procedure multiplies a 3x3 matrix and a 3x1 vector together.
        *
        *  author        : david vallado                  719-573-2600    1 mar 2001
        *
        *  inputs          description                    range / units
        *    mat         - 3 x 3 matrix
        *    vec         - vector
        *
        *  outputs       :
        *    vecout      - vector result of mat * vec
        *
        *  locals        :
        *    row         - row index
        *    col         - column index
        *    ktr         - index
        *
        *  coupling      :
        * --------------------------------------------------------------------------- */

            void    matvecmult
                (
                //std::vector< std::vector<double> > mat,
                          double mat[3][3],
                double vec[3],
                double vecout[3]
                )
            {
                int row, ktr;

                for (row = 0; row <= 2; row++)
                {
                    vecout[row] = 0.0;
                    for (ktr = 0; ktr <= 2; ktr++)
                        vecout[row] = vecout[row] + mat[row][ktr] * vec[ktr];
                }
            }  // matvecmult

    /* -----------------------------------------------------------------------------
    *
    *                           procedure angle
    *
    *  this procedure calculates the angle between two vectors.  the output is
    *    set to 999999.1 to indicate an undefined value.  be sure to check for
    *    this at the output phase.
    *
    *  author        : david vallado                  719-573-2600    1 mar 2001
    *
    *  inputs          description                    range / units
    *    vec1        - vector number 1
    *    vec2        - vector number 2
    *
    *  outputs       :
    *    theta       - angle between the two vectors  -Pi to Pi
    *
    *  locals        :
    *    temp        - temporary real variable
    *    magv1       - magnitude of vec1
    *    magv2       - magnitude of vec2
    *    small       - value defining a small value
    *    undefined   - large number to use in place of a not defined number
    *
    *  coupling      :
    *    dot           dot product of two vectors
    *    acos          arc cosine function
    *    mag           magnitude of a vector
    * --------------------------------------------------------------------------- */

        double  angle
            (
            double vec1[3],
            double vec2[3]
            )
        {
            double small, magv1, magv2, temp;
            small = 0.00000001;

            magv1 = mag(vec1);
            magv2 = mag(vec2);

            if (magv1*magv2 > small*small)
            {
                temp = dot(vec1, vec2) / (magv1*magv2);
                if (fabs(temp) > 1.0)
                    temp = sgn(temp) * 1.0;
                return acos(temp);
            }
            else
                return undefined;
        }  // angle
/* --------------------------------------------------------------------------- */
void sph2cart(double llh[3], double xyz[3])
{
    double cl=cos(llh[0]);
    double sl=sin(llh[0]);
    double cb=cos(llh[1]);
    double sb=sin(llh[1]);
    xyz[0]=llh[2]*cl*cb;
    xyz[1]=llh[2]*sl*cb;
    xyz[2]=llh[2]*sb;
}
void mattrans(double mat[3][3],double matt[3][3])
{
    for(int i = 0;i<3;i++)
        for(int j=0;j<3;j++)
            matt[i][j] = mat[j][i];
}
void cart2sph(double xyz[3], double llh[3])
{
    llh[2] = mag(xyz);
    llh[0] = atan2(xyz[1],xyz[0]);
    llh[1] = asin(xyz[2]/llh[2]);
}
/* ---------------------------------------------------------------------------
 * satbeam
 * input: satlon(deg), satlat(deg), sath(km), -- satellite position
 *        tarlon(deg), tarlat(deg), -- target position
 *        alf(deg) -- beam half angle
 * output: if satellite cover target, beam point( lon(deg), lat(deg) )
 *         else msg
--------------------------------------------------------------------------- */
#define D2R 0.0174532925199433
#define R2D 57.2957795130823
vector<vector<double> > SatBeam(vector<double>spos, double alf, double tlon, double tlat,double elev)
{
    vector<vector<double> > cov;
    //alf = alf*D2R;
    double x[3]={1,0,0};
    double z[3]={0,0,1};
    double rs = spos[2]/6378140;//(spos[2]+6378.140)/6378.140;
    double rc = sqrt(rs*rs-1);
    double llh[3]={spos[0]*D2R,spos[1]*D2R,rs};
    double xyz[3];
    sph2cart(llh,xyz);
    //printf("s: %f %f %f -> %f %f %f\n",llh[0]*180/pi,llh[1]*180/pi,llh[2],xyz[0],xyz[1],xyz[2]);
    double llht[3]={tlon*D2R,tlat*D2R,1};
    double xyzt[3];
    sph2cart(llht,xyzt);
    //printf("t: %f %f %f -> %f %f %f\n",llht[0]*180/pi,llht[1]*180/pi,llht[2],xyzt[0],xyzt[1],xyzt[2]);
    //printf("mag(s) = %f\n",mag(xyz));
    //printf("mag(t) = %f\n",mag(xyzt));
    double alfts = angle(xyz,xyzt);
    double beta = atan2(rs*cos(alfts)-1,rs*sin(alfts));
    //if(alfts > acos(1/rs))
    if( beta < elev*D2R)
    {
        //printf("can't cov target\n");
        return cov;
    }
    double alfzs= (90 -spos[1])*D2R;//angle(xyz,z);
    //printf("alfts = %f alfzs = %f\n",alfts*180/pi,alfzs*180/pi);
    //double cx[3],cy[3],cz[3];
    double matse[3][3],mates[3][3],mates1[3][3],matse1[3][3];
    double ts[3];
    addvec(-1,xyz,1,xyzt,ts);
    //double rts = mag(ts);
    //printf("rts = %f\n",rts);
    norm(ts,mates[2]);//norm(ts,cz);
    if(alfts > 5*D2R)
        cross(mates[2],xyz,mates[1]);//cross(cz,xyz,cy);
    else if(alfzs > 5*D2R)
        cross(z,xyz,mates[1]);//cross(z,xyz,cy);
    else
        cross(x,xyz,mates[1]);//cross(x,xyz,cy);
    norm(mates[1],mates1[1]);//norm(cy,cy);
    memcpy(mates[1],mates1[1],sizeof(double)*3);
    cross(mates[1],mates[2],mates[0]);//cross(cy,cz,cx);
    mattrans(mates,matse);
    //printf("Ces:\n%f %f %f\n%f %f %f\n%f %f %f\n",mates[0][0],mates[0][1],mates[0][2],
    //        mates[1][0],mates[1][1],mates[1][2],mates[2][0],mates[2][1],mates[2][2]);

    norm(xyz,mates1[2]);
    cross(mates1[1],mates1[2],mates1[0]);
    mattrans(mates1,matse1);
    //printf("%f %f %f\n%f %f %f\n%f %f %f\n",cx[0],cx[1],cx[2],cy[0],cy[1],cy[2],cz[0],cz[1],cz[2]);
    //printf("Ces1:\n%f %f %f\n%f %f %f\n%f %f %f\n",mates1[0][0],mates1[0][1],mates1[0][2],
    //        mates1[1][0],mates1[1][1],mates1[1][2],mates1[2][0],mates1[2][1],mates1[2][2]);
    double so[3];
    matvecmult(mates,xyz,so);
    //printf("so = %f %f %f\n",so[0],so[1],so[2]);
    //double gama = atan2(so[0],-so[2]);
    //printf("gama = %f\n",gama*180/pi);
    double azc, stp;
    double ccl[3];
    ccl[0] = (rc+cos(alf)*so[2])/(sin(alf)*so[0]);
    if(ccl[0] > -1)
    {
        ccl[0] = acos(-ccl[0]);
        ccl[1] = pi/2 - alf;
        ccl[2] = rc;
        double ccls[3], ccle[3], ccls1[3];
        sph2cart(ccl,ccls);
        matvecmult(matse,ccls,ccle);
        matvecmult(mates1,ccle,ccls1);
        azc = atan2(ccls1[1],-ccls1[0]);
        int n = ccl[0]/(5*D2R);
        stp = azc/n;
        //printf("azc = %f %f %f\n",azc*180/pi,ccl[0]*180/pi,stp*180/pi);
    }
    double alfc = asin(1/rs);
    double cms[3],cmc[3],cmce[3],xyzc[3],llhc[3];
    int j = 0;
    for(int i=0;i<72;i++)
    {
        cms[0] = (i*5-180)*D2R;
        double b = 2*(cos(alf)*so[2]+sin(alf)*so[0]*cos(cms[0]));
        double dlt = b*b - 4*(rs*rs-1);
        //printf("b=%f,dlt=%f\n",b,dlt);
        //double rp;
        if(dlt>0)
        {
            cms[1] = pi/2-alf;
            cms[2]=-(b+sqrt(dlt))/2;
            sph2cart(cms,cmc);
            matvecmult(matse,cmc,cmce);
        }
        else
        {
            cms[0] = pi+azc - j*stp;
            if(cms[0] < pi-azc)
                cms[0] = pi-azc;
            cms[1] = -pi/2 + alfc;
            cms[2] = rc;
            //printf("%d %d %f %f %f\n",i,j,cms[0]*180/pi,cms[1]*180/pi,cms[2]);
            sph2cart(cms,cmc);
            //printf("%d %d %f %f %f\n",i,j,cmc[0],cmc[1],cmc[2]);
            matvecmult(matse1,cmc,cmce);
            j++;
            /*cms[2]=rc;
            double a = k+cos(gama);
            b=2*cos(cms[0])*sin(gama);
            double c=k-cos(gama);
            dlt = b*b - 4*a*c;
            cms[1] = pi/2 - 2*atan2(-b+sqrt(dlt),2*a);*/

        }
        //sph2cart(cms,cmc);
        //matvecmult(matse,cmc,cmce);
        //printf("%f %f %f\n",cmce[0],cmce[1],cmce[2]);
        addvec(1,cmce,1,xyz,xyzc);
        cart2sph(xyzc,llhc);
        //printf("%f,%f,\n",llhc[0]*R2D,llhc[1]*R2D);
        vector<double>ll;
        ll.push_back(llhc[0]*R2D);
        ll.push_back(llhc[1]*R2D);
        cov.push_back(ll);
    }
    return cov;
}
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
// 创建可旋转的帆板
osg::ref_ptr<osg::MatrixTransform> createRotatingSail(double angle) {
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;

    // 帆板顶点数据
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->push_back(osg::Vec3(0, -1, 0.0f));
    vertices->push_back(osg::Vec3(-1, -2, 0.0f));
    vertices->push_back(osg::Vec3(-1, -5, 0.0f));
    vertices->push_back(osg::Vec3(1, -5, 0.0f));
    vertices->push_back(osg::Vec3(1, -2, 0.0f));
    geometry->setVertexArray(vertices);
    // 帆板颜色
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(0.5f, 0.5f, 0.5f, 1.0f));  // 白色
    geometry->setColorArray(colors);
    geometry->setColorBinding(osg::Geometry::BIND_OVERALL);

    // 帆板索引
    osg::ref_ptr<osg::DrawElementsUShort> indices = new osg::DrawElementsUShort(osg::PrimitiveSet::POLYGON);
    for (int i = 0; i < 5;++i){//vertices->getNumElements()+2; ++i) {
        indices->push_back(i);
    }
    geometry->addPrimitiveSet(indices);

    osg::ref_ptr<osg::Geometry> geometry1 = new osg::Geometry;
    osg::ref_ptr<osg::Vec3Array> vertices1 = new osg::Vec3Array;
    vertices1->push_back(osg::Vec3(0,  1, 0.0f));
    vertices1->push_back(osg::Vec3(1,  2, 0.0f));
    vertices1->push_back(osg::Vec3(1,  5, 0.0f));
    vertices1->push_back(osg::Vec3(-1, 5, 0.0f));
    vertices1->push_back(osg::Vec3(-1, 2, 0.0f));
    geometry1->setVertexArray(vertices1);
    geometry1->setColorArray(colors);
    geometry1->setColorBinding(osg::Geometry::BIND_OVERALL);
    geometry1->addPrimitiveSet(indices);

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geometry);
    geode->addDrawable(geometry1);

    osg::ref_ptr<osg::MatrixTransform> transform = new osg::MatrixTransform;
    transform->setMatrix(osg::Matrix::rotate(angle, 0, 1, 0));
    transform->addChild(geode);

    /*/ 旋转回调类
    class RotateCallback : public osg::NodeCallback {
    public:
        RotateCallback() : angle(0.0) {}
        virtual void operator()(osg::Node* node, osg::NodeVisitor* nv) {
            osg::MatrixTransform* mt = dynamic_cast<osg::MatrixTransform*>(node);
            if (mt) {
                angle += 0.001;
                mt->setMatrix(osg::Matrix::rotate(angle, 0, 1, 0));
            }
            traverse(node, nv);
        }
    private:
        double angle;
    };

    transform->setUpdateCallback(new RotateCallback);*/
    return transform;
}
// 创建六面颜色各异的六面体
osg::ref_ptr<osg::Node> createColoredCube() {
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;

    // 顶点数据
    osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
    vertices->push_back(osg::Vec3(-1.0f, -1.0f, -1.0f));
    vertices->push_back(osg::Vec3(1.0f, -1.0f, -1.0f));
    vertices->push_back(osg::Vec3(1.0f, 1.0f, -1.0f));
    vertices->push_back(osg::Vec3(-1.0f, 1.0f, -1.0f));
    vertices->push_back(osg::Vec3(-1.0f, -1.0f, 1.0f));
    vertices->push_back(osg::Vec3(1.0f, -1.0f, 1.0f));
    vertices->push_back(osg::Vec3(1.0f, 1.0f, 1.0f));
    vertices->push_back(osg::Vec3(-1.0f, 1.0f, 1.0f));
    geometry->setVertexArray(vertices);

    // 颜色数据，每个面一种颜色
    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    colors->push_back(osg::Vec4(1.0f, 0.0f, 0.0f, 1.0f));  // 红色
    colors->push_back(osg::Vec4(0.0f, 1.0f, 0.0f, 1.0f));  // 绿色
    colors->push_back(osg::Vec4(0.0f, 0.0f, 1.0f, 1.0f));  // 蓝色
    colors->push_back(osg::Vec4(1.0f, 1.0f, 0.0f, 1.0f));  // 黄色
    colors->push_back(osg::Vec4(1.0f, 0.0f, 1.0f, 1.0f));  // 紫色
    colors->push_back(osg::Vec4(0.0f, 1.0f, 1.0f, 1.0f));  // 青色
    geometry->setColorArray(colors);
    geometry->setColorBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);

    // 索引数据，定义每个面的顶点
    const unsigned int faceIndices[6][4] = {
        {0, 1, 2, 3},  // 前面
        {4, 5, 6, 7},  // 后面
        {0, 3, 7, 4},  // 左面
        {1, 2, 6, 5},  // 右面
        {0, 1, 5, 4},  // 底面
        {2, 3, 7, 6}   // 顶面
    };
    for (int i = 0; i < 6; ++i) {
        osg::ref_ptr<osg::DrawElementsUShort> indices = new osg::DrawElementsUShort(osg::PrimitiveSet::QUADS);
        for (int j = 0; j < 4; ++j) {
            indices->push_back(faceIndices[i][j]);
        }
        geometry->addPrimitiveSet(indices);
    }

    osg::ref_ptr<osg::Geode> geode = new osg::Geode;
    geode->addDrawable(geometry);
    return geode;
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
/*
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
osg::Vec3d SatEci(const cSatellite& sat, time_t t)
{
    cJulian cjt(t);
    cEciTime sateci = sat.PositionEci(cjt);
    cEci pos = sateci;
    osg::Vec3d xyz(pos.Position().m_x,pos.Position().m_y,pos.Position().m_z);
    return xyz;
}
// Reference time for the J2000 ECI coordinate frame
static DateTime J2000Epoch(2000, 1, 1, 12.00);

// Transform that takes us from a J2000 ECI reference frame
// to an ECEF reference frame (i.e. MapNode)
class J2000ToECEFTransform : public osg::MatrixTransform
{
public:
    void setDateTime(const DateTime& dt)
    {
        osg::Matrix matrix = createMatrix(dt);
        setMatrix(matrix);
    }

    static osg::Matrix createMatrix(const DateTime& dt)
    {
        // Earth's rotation rate: International Astronomical Union (IAU) GRS 67
        const double IAU_EARTH_ANGULAR_VELOCITY = 7292115.1467e-11; // (rad/sec)

        double secondsElapsed = (double)(dt.asTimeStamp() - J2000Epoch.asTimeStamp());
        const double rotation = IAU_EARTH_ANGULAR_VELOCITY * secondsElapsed;

        osg::Matrix matrix;
        matrix.makeRotate(rotation, 0, 0, 1);
        return matrix;
    }
};
//    osg::Matrix eci2ecef = J2000ToECEFTransform::createMatrix(t0);
//        osg::Vec3d eci = SatEci(*sat,t.asTimeStamp());
//        printf("eci: %lf %lf %lf\n",eci.x(),eci.y(),eci.z());
//	osg::Vec3d ecef =eci2ecef*eci;
//        printf("ecef: %lf %lf %lf\n",ecef.x(),ecef.y(),ecef.z());
*/
