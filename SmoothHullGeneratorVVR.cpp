#include "SmoothHullGeneratorVVR.h"
#include <iostream>
#include <fstream>
#include <limits>
//#include <process.h>


//#define DISPLAY_INFO

namespace SCD
{


  SmoothHullGeneratorVVR::SmoothHullGeneratorVVR(double r, double R):
    _epsilon(1e-8)
  {
    _r = r;
    _R = R;
  }


  bool SmoothHullGeneratorVVR::findCenter(int p1, int p2, int p3, vector3d &center)
  {
    /**
     *we will solve the following system :
     *	(1)	u.x = u.u/2			(median plan of [p1p2])
     *	(2)	v.x = v.v/2			(median plan of [p1p3])
     *	(3)	x.x = (R-r)^2		(norm of x)
     *	(4)	x.n < 0				(discriminant between the two solution)
     *where u is the vector p1p2, v is p1p3, n is the normal to p1p2p3 and x is the unknown vector p1center
     *System is solved writing x in the (u,v,u*v) frame : x = au + bv + c(u*v)
     *(1) becomes au.u+bu.v = u.u/2 and (2) au.v+bv.v = v.v/2. Both equations are solved together
     *c is then compute from (3) while respecting (4) which is equivalent to x.(u*v)<0 <=> c<0
     */
    double a,b,c,duv,nu,nv,nuv;
    vector3d u,v,uv;
    u = _points[p2]-_points[p1];
    v = _points[p3]-_points[p1];
    duv = u%v;		//dot product
    uv = u^v;		//cross product
    nu = u.normsquared();
    nv = v.normsquared();
    nuv = uv.normsquared();

    if (nuv == 0) { // aligned points
#ifdef DISPLAY_INFO
      std::cout << "WARNING, no sphere can fit on face made of 3 aligned points" << std::endl;
#endif
      return false;
    }
    a = (nu*nv-nv*duv)/(2*nuv);
    b = (nu*nv-nu*duv)/(2*nuv);
    if (((_R-_r)*(_R-_r)-a*a*nu-b*b*nv-2*a*b*duv)/nuv<=0) {
#ifdef DISPLAY_INFO
      //std::cout << "WARNING, no sphere of radius " << _R << " can fit on face "  << p1 << " " << p2 << " " << p3 << std::endl;
#endif
      return false;
    } else {
      c = -sqrt(((_R-_r)*(_R-_r)-a*a*nu-b*b*nv-2*a*b*duv)/nuv);

      center = u*a + v*b + uv*c + _points[p1];
      return true;
    }
  }


  bool SmoothHullGeneratorVVR::isInSphere(vector3d &point, vector3d &center)
  {
    vector3d v(point,center);
    bool res = (v.normsquared() <= ((_R - _r) * (_R - _r) * (1.+_epsilon)));
    return res;
  }

  double SmoothHullGeneratorVVR::distInSphere(vector3d &point, vector3d &center)
  {
    vector3d v(point,center);
    return (v.normsquared() - ((_R - _r) * (_R - _r)));// * (1.+_epsilon)));
  }
  bool SmoothHullGeneratorVVR::allPointsInSphere(vector3d &center)
  {
    bool b = true;
    for (unsigned int l=0; l<_points.size() && b; ++l)
    {
      b = isInSphere(_points[l],center);
    }
    return b;
  }

  bool SmoothHullGeneratorVVR::findFirstTriangle(unsigned &i,unsigned &j,unsigned &k, vector3d &c)
  {
    unsigned int n = _points.size();

    bool b = false;
    for (i=0; i<n; ++i)
    {
      for (j=i+1; j<n; ++j)
      {
        for (k=j+1; k<n; ++k)
        {
          //let's test if every point is inside the sphere [i,j,k]
          if (findCenter(i, j, k, c))
            b = allPointsInSphere(c);
          if (b) break;
          //now let's test if every point is inside the sphere [i,k,j]
          if (findCenter(i, k, j, c))
            b = allPointsInSphere(c);
          if (b)
          {
            _ccw = false;
            break;
          }
        }
        if (b) break;
      }
      if (b) break;
    }
    return b;
  }

  double SmoothHullGeneratorVVR::distMaxPointsInSphere(vector3d &center)
  {
    double err;
    double err_max = -(_R-_r)*(_R-_r);
    for (unsigned int l=0; l<_points.size(); ++l)
    {
      err = distInSphere(_points[l],center);
      if (err>err_max)
        err_max = err;
      if (err<=-(_R-_r)*(_R-_r)) {
        std::cout << "PROBLEM, point at center ? " << err << " center " << center << " point " << _points[l] << std::endl;
      }
    }
    if (err_max==-(_R-_r)*(_R-_r)) {
      std::cout << "PROBLEM, all points at center ? " << " center " << center << " nb points " << _points.size() << std::endl;
    }
    //std::cout << err_max << std::endl;
    return err_max;
  }

  double SmoothHullGeneratorVVR::computeAngle(vector3d &p1, vector3d &p2, vector3d &p3, vector3d &axe)
  {
    vector3d u(p2,p1);
    vector3d v(p2,p3);
    double cosinus = u%v/(u.norm()*v.norm());
    if (cosinus > double(1)) {
#ifdef DISPLAY_INFO
      std::cout << "ERROR in computing angle, found cos(angle)>1 " << cosinus << std::endl;
#endif
      cosinus = 1;
      //return 0;
    } else if (cosinus < double(-1)) {
#ifdef DISPLAY_INFO
      std::cout << "ERROR in computing angle, found cos(angle)<-1 " << cosinus << std::endl;
#endif
      cosinus = -1;
    }
    return (((axe%(u^v)>double(0))?1:-1)*acos(cosinus));
  }
  //
  void SmoothHullGeneratorVVR::cover(void)
  {
    unsigned int n = _points.size();


    //look for largest distance
    vector3d d;
    double d_max = 0;
    for (unsigned i=0; i<n; ++i)
    {
      for (unsigned j=i+1; j<n; ++j)
      {
        d = _points[i] - _points[j];
        if (d_max<d.norm())
          d_max = d.norm();
      }
    }
    std::cout << "Maximum body distance " << d_max << std::endl;
    if (d_max>=2.*(_R-_r))
      std::cout << "ERROR, impossible to compute an STP-BV for this body with R=" << _R << ", choose larger R" << std::endl;

    // 1 - find a triangle to start with
    std::cout << "Searching for initial triangle... ";
    double epsilon_temp = _epsilon;
    _epsilon = 0.; // should find an initial triangle with such epsilon
    bool b = false;
    _ccw = true;			//is the triangle describes by [i,j,k] counterclockwise
    unsigned int i,j,k;
    vector3d c;

#if 1  // find the first matching triangle
    if (!(b=findFirstTriangle(i,j,k,c)))
    {
      _epsilon = 1e-8;
      std::cout << "Failed" << std::endl <<  "trying with more compliant constraints... " ;
      b=findFirstTriangle(i,j,k,c);
    }

#else  // find the triangle with maximum distance of other vertices inside, better but can be very time consuming for large clouds
    int i_min = -1;
    int j_min = -1;
    int k_min = -1;
    double dist_min = 100000.;
    double dist;
    vector3d c_min;
    for (i=0; i<n; ++i)
    {
      for (j=i+1; j<n; ++j)
      {
        for (k=j+1; k<n; ++k)
        {
          //let's test if every point is inside the sphere [i,j,k]
          if (findCenter(i, j, k, c)) {
            dist = distMaxPointsInSphere(c);
            if (dist<dist_min) {
              dist_min = dist;
              i_min = i;
              j_min = j;
              k_min = k;
              c_min = c;
              _ccw = true;
            }
          }
          //now let's test if every point is inside the sphere [i,k,j]
          if (findCenter(i, k, j, c)) {
            dist = distMaxPointsInSphere(c);
            if (dist<dist_min) {
              dist_min = dist;
              i_min = i;
              j_min = j;
              k_min = k;
              c_min = c;
              _ccw = false;
            }
          }
        }
      }
      std::cout << i << "/" << n << std::endl;
    }
    i = i_min;
    j = j_min;
    k = k_min;
    c = c_min;
    std::cout << "best initial triangle " << i << " " << j << " " << k << " dist " << dist_min << std::endl;
#endif
    _epsilon = epsilon_temp;

    if (!b)		//no initial triangle found
    {
      std::cout << "failed" << std::endl;
      return;
    }

    if (!_ccw) 
    {
      unsigned int e = j;
      j = k;
      k = e;
    }
    std::cout << "found" << std::endl;
    //		std::cout << i << ", " << j << ", " << k << std::endl;
    VVRSFace s = {i,j,k,c};
    _index.insert(i);
    _index.insert(j);
    _index.insert(k);
    //		printSphere(s);

    // 2 - turn around the edge
    std::cout << "Computing hull... ";
    travelCover(s);
    if (_spheres.size() != ((_index.size()-2)*2))
    {
      std::cout << "WARNING : EULER FORMULA IS NOT RESPECTED" << std::endl;
    }
    else
    {
      std::cout << "done" << std::endl;
    }
  }

  int	SmoothHullGeneratorVVR::getKey(int a, int b)
  {
    return ((a<b)?(a*_points.size()+b):(b*_points.size()+a));
  }

  int	SmoothHullGeneratorVVR::getKey(turnData& td)
  {
    return ((td._p1<td._p2)?(td._p1*_points.size()+td._p2):(td._p2*_points.size()+td._p1));
  }

  //
  double SmoothHullGeneratorVVR::getKeyByAngle(turnData& td)
  {
    double aMin = 100.;
    for (unsigned int j=0; j<_points.size(); ++j)
    {
      if (j!=td._p1 && j!=td._p2 && j!=td._p3)
      {
        vector3d c;
        findCenter(td._p1, j, td._p2, c);
        if (allPointsInSphere(c))
        {
          vector3d pointp1p2moy((_points[td._p1]+_points[td._p2])/2);
          vector3d pointp1mp2(_points[td._p2]-_points[td._p1]);
          double a = fabs(computeAngle(td._previousCenter,pointp1p2moy,c,pointp1mp2));
          if (a<aMin)
            aMin=a;
        }
      }
    }
    return aMin;
  }

  void SmoothHullGeneratorVVR::travelCover(VVRSFace s)
  {
    std::set<int> edgeStack;
    std::multimap<double,turnData> edgeStackByAngle;
    std::set<int> computedEdge;
    turnData td = {s._point1,s._point2,s._point3,s._center};
    edgeStack.insert(getKey(td));
    edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td),td));
    td._p1 = s._point2;
    td._p2 = s._point3;
    td._p3 = s._point1;
    edgeStack.insert(getKey(td));
    edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td),td));
    td._p1 = s._point3;
    td._p2 = s._point1;
    td._p3 = s._point2;
    edgeStack.insert(getKey(td));
    edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td),td));

    _spheres.push_back(s);
    _index.insert(s._point1);
    _index.insert(s._point2);
    _index.insert(s._point3);

#ifdef DISPLAY_INFO
    std::cout.precision(17);
    std::cout << "Original triangle : " << s._point1 << " " << s._point2 << " " << s._point3 << std::endl;
    std::cout << "Center " << s._center << std::endl;
    std::cout << _points[s._point1] << std::endl;
    std::cout << _points[s._point2] << std::endl;
    std::cout << _points[s._point3] << std::endl;
    for (std::multimap<double, turnData>::iterator ei=edgeStackByAngle.begin(); ei!=edgeStackByAngle.end(); ++ei)
      std::cout << "Edge [" << (ei->second)._p1 << "," << (ei->second)._p2 << "] angle " << ei->first << std::endl;
#endif
    double curr_edge_angle;

    while (!edgeStackByAngle.empty())
    {
postponed:
      td = (edgeStackByAngle.begin())->second;
      //#ifdef DISPLAY_INFO
      curr_edge_angle = (edgeStackByAngle.begin())->first;
      //#endif
      edgeStack.erase(getKey(td));
      edgeStackByAngle.erase(edgeStackByAngle.begin());

      if (!((computedEdge.insert(getKey(td))).second))
      {
        std::cout << "WARNING, edge already processed : "<< td._p1 << ", " << td._p2 << std::endl;
        continue;
      }



      std::multimap<double, int> sortedPoints;
      vector3d c;
      unsigned int j_min;
      double distInSphereMin = 100000.;
      for (unsigned int j=0; j<_points.size(); ++j)
      {
        if (j!=td._p1 && j!=td._p2 && j!=td._p3)
        {
          if (findCenter(td._p1, j, td._p2, c)) {
            double err_max = distMaxPointsInSphere(c);
            if (err_max<distInSphereMin) {
              distInSphereMin = err_max;
              j_min = j;
            }
#ifdef DISPLAY_INFO
            std::cout << "err points in sphere " << j << " " << err_max << std::endl;
#endif
            if (allPointsInSphere(c))
            {
              vector3d pointp1p2moy((_points[td._p1]+_points[td._p2])/2);
              vector3d pointp1mp2(_points[td._p2]-_points[td._p1]);
              double a = computeAngle(td._previousCenter,pointp1p2moy,c,pointp1mp2);
              sortedPoints.insert(std::pair<double, int>(a,j));
            }
          }
        }
      }
      if (sortedPoints.empty()) {
        std::cout << "WARNING, no new vertex found in rotation around edge [" << td._p1 << ", " << td._p2 << "] that satisfy inclusion in shpere of all other vertices with chosen precision" << std::endl;
        std::cout << "Taking vertex that satisfy best this condition with value " << distInSphereMin << std::endl;
        if (!findCenter(td._p1, j_min, td._p2, c))
          std::cout << "ERROR, choosen new vertex does not allow to build a shpere with edge" << std::endl;
        vector3d pointp1p2moy((_points[td._p1]+_points[td._p2])/2);
        vector3d pointp1mp2(_points[td._p2]-_points[td._p1]);
        double a = computeAngle(td._previousCenter,pointp1p2moy,c,pointp1mp2);
        sortedPoints.insert(std::pair<double, int>(a,j_min));
      }
      if (!sortedPoints.empty())
      {
        int p = sortedPoints.begin()->second;
#ifdef DISPLAY_INFO
        std::cout << "nb candidate vertices " << sortedPoints.size() << " around edge [" << td._p1 << ", " << td._p2 << "] ";
        for (std::multimap<double, int>::iterator pi=sortedPoints.begin(); pi!=sortedPoints.end(); ++pi) {
          std::cout << " " << pi->first << " " << pi->second << " " << (_index.find(pi->second)!=_index.end()); //->first;
          if (pi->second>1.) {
            std::cout << " computed_edge_angle " << curr_edge_angle << " next " <<  (edgeStackByAngle.begin())->first; //->first;
          }
        }
        std::cout << std::endl;
        // see if none of the possible new edges already computed
        if (computedEdge.find(getKey(p,td._p1))!=computedEdge.end())
          std::cout << "EDGE with " << p << " and " << td._p1 << " already added" << std::endl;
        if (computedEdge.find(getKey(p,td._p2))!=computedEdge.end())
          std::cout << "EDGE with " << p << " and " << td._p2 << " already added" << std::endl;
#endif
        if (sortedPoints.size()>1) {
#ifdef DISPLAY_INFO
          std::cout << td._p1 << " " << _points[td._p1] << std::endl << td._p2 << " " << _points[td._p2] << std::endl;
          for (std::multimap<double, int>::iterator pi=sortedPoints.begin(); pi!=sortedPoints.end(); ++pi)
            std::cout << pi->second << " " << _points[pi->second] << std::endl; //->first;
#endif
          // if there is indeterminacy between candidate, postpone choice by increasing angle
          if (curr_edge_angle<3.14159) {
            computedEdge.erase(getKey(td));
            edgeStack.insert(getKey(td));
            edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td)+3.14159,td));
#ifdef DISPLAY_INFO
            std::cout << "POSTPONING decision" << std::endl;
#endif
            goto postponed;
#ifdef DISPLAY_INFO
          } else {
            std::cout << "treating POSTPONED edge" << std::endl;
#endif
          }
        }
        vector3d cs;
        findCenter(td._p1, p, td._p2, cs);
        VVRSFace sp = {td._p1, p, td._p2,cs};
        _spheres.push_back(sp);
        _index.insert(p);


        turnData td1 = {td._p1,p, td._p2,cs};
        turnData td2 = {p, td._p2,td._p1,cs};


        if (!(edgeStack.insert(getKey(td1)).second))
        {
          std::multimap<double,turnData>::iterator it = edgeStackByAngle.begin();
#ifdef DISPLAY_INFO
          std::cout << "erasing edge "<< td1._p1 << " " << td1._p2 << " completely processed" << std::endl;
#endif
          while (getKey(it->second) != getKey(td1))
            ++it;
          edgeStackByAngle.erase(it);
          edgeStack.erase(getKey(td1));
        }
        else
        {
          edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td1),td1));
        }

        if (!(edgeStack.insert(getKey(td2)).second))
        {
          std::multimap<double,turnData>::iterator it = edgeStackByAngle.begin();
#ifdef DISPLAY_INFO
          std::cout << "erasing edge "<< td2._p1 << " " << td2._p2 << " completely processed" << std::endl;
#endif
          while (getKey(it->second) != getKey(td2))
            ++it;
          edgeStackByAngle.erase(it);
          edgeStack.erase(getKey(td2));
        }
        else
        {
          edgeStackByAngle.insert(std::pair<double,turnData>(getKeyByAngle(td2),td2));
        }

      }
      else
      {
        std::cout << "No sphere were found while turning around [" << td._p1 << ", " << td._p2 << "]" << std::endl;
      }
    }

#ifdef DISPLAY_INFO
    for (int i=0; i<_spheres.size(); ++i) {
      std::cout << _spheres[i]._point1 << " " << _spheres[i]._point2 << " " << _spheres[i]._point3 << std::endl;
    }
    for (unsigned int j=0; j<_points.size(); ++j) {
      std::cout << j << " " << _points[j].x << " " << _points[j].y << " " << _points[j].z << std::endl;
    }
#endif
  }

  //#ifdef FALSE
  void SmoothHullGeneratorVVR::printSphere(VVRSFace& s)
  {
    std::cout << "**Spheres**" << std::endl;
    std::cout<< "   center : " << s._center.x << ", " << s._center.y << ", " << s._center.z << std::endl;
    std::cout<< "   point1 [" << s._point1 << "] : " << _points[s._point1].x << ", " << _points[s._point1].y << ", " << _points[s._point1].z << std::endl;
    std::cout<< "   point2 [" << s._point2 << "] : " << _points[s._point2].x << ", " << _points[s._point2].y << ", " << _points[s._point2].z << std::endl;
    std::cout<< "   point3 [" << s._point3 << "] : " << _points[s._point3].x << ", " << _points[s._point3].y << ", " << _points[s._point3].z << std::endl;
  }

  void SmoothHullGeneratorVVR::readVertex(const std::string& filename)
  {
    FILE * file;
    // Open file
    file = fopen(filename.c_str(), "r");
    if (file == 0)
    {
      std::cout << "unable to open file " << filename << std::endl;
      return;
    }

    std::cout << "reading " << filename << "..... ";
    int nbPtsPerFace;
    if(fscanf(file, "%d\n", &nbPtsPerFace) != 1)
    {
      std::cout << filename << " is in the wrong file format." << std::endl;
      return;
    }

    int nVert;
    double x,y,z;
    _points.clear();
    if(fscanf(file, "%d\n", &nVert) != 1)
    {
      std::cout << filename << " is in the wrong file format." << std::endl;
      return;
    }

    vector3d v;
    for (int i=0; i<nVert; ++i)
    {
      if(fscanf(file, "%lf %lf %lf\n", &x, &y, &z) != 3)
      {
        std::cout << filename << " is in the wrong file format." << std::endl;
        return;
      }

      v.x = double(x);
      v.y = double(y);
      v.z = double(z);
      _points.push_back(v);
    }
  }

  void SmoothHullGeneratorVVR::output(const std::string& rootPath)
  {
    //we first build a map to change the index of points in the spheres, according to the vertex that have been kept
    int k=0;
    std::vector<int> map;

    std::cout << "please enter the name of the output file" << std::endl;
    std::string s;
    std::cin >> s;
    std::string name = rootPath + s;

    std::ofstream output;
    output.open(name.c_str(),std::ofstream::out);
    if (output.is_open())
    {
      output << _r << "," << _R << "\n";
      output << _index.size() << "," << _spheres.size() << "\n";
      for (unsigned int i=0; 0<_index.size();++i)
      {
        if (*(_index.begin()) == i)
        {
          output << "[" << _points[*(_index.begin())].x << ","
            << _points[*(_index.begin())].y << ","
            << _points[*(_index.begin())].z << "],";
          map.push_back(k);
          _index.erase(_index.begin());
          k++;
        }
        else
        {
          map.push_back(-1);
        }
      }
      output << std::endl;
      for (unsigned int i=0; i<_spheres.size(); ++i)
      {
        output << "[" << (map[_spheres[i]._point1]+1) << "," << (map[_spheres[i]._point2]+1) << "," << (map[_spheres[i]._point3]+1) << "],";
      }
      output.close();
    }
  }
  //#endif

  void SmoothHullGeneratorVVR::loadGeometry(const std::string& filename)
  {
    readVertex(filename);
  }
  //#ifdef FALSE

  void SmoothHullGeneratorVVR::compute3DSMaxHull(const std::string& rootPath)
  {
    _index.clear();
    _spheres.clear();
    cover();
    output(rootPath);
  }

  void SmoothHullGeneratorVVR::computeSmoothHull(std::vector<vector3d>& outPoints, std::vector<VVRSFace>& outSFaces)
  {
    //computing hull
    _index.clear();
    _spheres.clear();
    cover();

    //output result in the return parameters
    std::set<int> computedEdge;
    outPoints.clear();
    outPoints.reserve(_index.size());
    outSFaces.clear();
    outSFaces.reserve(_spheres.size());
    //		outEdges.clear();
    //		outEdges.resize(_index.size() + _spheres.size() -2);		//Euler Formula V-E+F = 2
    //		adjacentFaces.clear();
    //		adjacentFaces.reserve(_index.size());
    int k=0;
    std::vector<int> map;

    std::cout << "size of _index : " << _index.size() << std::endl;
    for (unsigned int i=0; 0<_index.size();++i)
    {
      if (*(_index.begin()) == i)
      {
        vector3d v;
        //				std::vector<int> vi;
        //				vi.clear();
        v = _points[*(_index.begin())];
        map.push_back(k);
        outPoints.push_back(v);
        //				adjacentFaces.push_back(vi);
        _index.erase(_index.begin());
        k++;
      }
      else
      {
        map.push_back(-1);
      }
    }
    for (unsigned int i=0; i<_spheres.size(); ++i)
    {
      VVRSFace sf;
      sf._point1 = map[_spheres[i]._point1];
      sf._point2 = map[_spheres[i]._point2];
      sf._point3 = map[_spheres[i]._point3];
      sf._center = _spheres[i]._center;
      outSFaces.push_back(sf);
    }
  }
  //#endif
  void SmoothHullGeneratorVVR::computeVVR_Prime(const std::string& filename)
  {
    std::ofstream os;
    std::vector<VVRsphere> smallSpheres;
    std::vector<VVRsphere> bigSpheres;
    std::vector< std::pair<int, VVRtorus> > torus;
    std::map<int, int> tinds;
    std::vector<std::vector<VVRcone> > ssVVR;
    std::vector<std::vector<bool> > computed(_points.size());
    std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > > tVVR;
    std::vector<SCD::faceVVR> bsVVR;
    VVRsphere stmp;
    VVRcone ctmp;
    double l;
    VVRtorus ttmp;
    int scount = 0;
    int tind, i;
    std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > tpair;
    SCD::faceVVR ftmp;
    double epsilon = 1e-10;

    //computing hull
    _index.clear();
    _spheres.clear();
    cover();

    std::vector<int> ind(_index.size());
    std::vector<int> invind(_points.size(),-1);

    int difference=_points.size()-_index.size();

    i=0;
    for (std::set<int,std::less<int> >::iterator it=_index.begin(); it!=_index.end() ; ++it,++i)
    {
      ind[i]=*it;
      invind[*it]=i;
    }

    for(unsigned(i) = 0 ; i < _points.size() ; ++i)
    {
      computed[i].resize(_points.size());
      for(unsigned j = 0 ; j < _points.size() ; ++j)
        computed[i][j] = false;
    }

    for(std::vector<vector3d>::const_iterator it = _points.begin() ; 
        it != _points.end() ; 
        ++it)
    {
      stmp._center = *it;
      stmp._radius = _r;
      smallSpheres.push_back(stmp);
    }
    ssVVR.resize(_points.size());

    for(std::vector<VVRSFace>::const_iterator it = _spheres.begin() ; 
        it != _spheres.end() ; 
        ++it)
    {
      //register the new big sphere
      stmp._center = (*it)._center;
      stmp._radius = _R;
      bigSpheres.push_back(stmp);

      tind = getKey((*it)._point1, (*it)._point2);
      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point1] + _points[(*it)._point2]) / 2;
        l = (_points[(*it)._point1] - _points[(*it)._point2]).norm();
        ttmp._normal = (_points[(*it)._point1] - _points[(*it)._point2]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane1.first = tinds[tind];
      ftmp._plane1.second = (_points[(*it)._point1] - stmp._center)^(_points[(*it)._point2] - stmp._center);
      ftmp._plane1.second /= ftmp._plane1.second.norm();

      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the VVR plane for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane1;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.first = ctmp;
        // [point2, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.second = ctmp;
        computed[(*it)._point1][(*it)._point2] = computed[(*it)._point2][(*it)._point1] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR plane for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane1;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point2, (*it)._point3);
      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point2] + _points[(*it)._point3]) / 2;
        l = (_points[(*it)._point2] - _points[(*it)._point3]).norm();
        ttmp._normal = (_points[(*it)._point2] - _points[(*it)._point3]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane2.first = tinds[tind];
      ftmp._plane2.second = (_points[(*it)._point2] - stmp._center)^(_points[(*it)._point3] - stmp._center);
      ftmp._plane2.second /= ftmp._plane2.second.norm();

      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane2;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point2, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.first = ctmp;
        // [point3, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.second = ctmp;
        computed[(*it)._point2][(*it)._point3] = computed[(*it)._point3][(*it)._point2] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane2;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point3, (*it)._point1);
      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point3] + _points[(*it)._point1]) / 2;
        l = (_points[(*it)._point3] - _points[(*it)._point1]).norm();
        ttmp._normal = (_points[(*it)._point3] - _points[(*it)._point1]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane3.first = tinds[tind];
      ftmp._plane3.second = (_points[(*it)._point3] - stmp._center)^(_points[(*it)._point1] - stmp._center);
      ftmp._plane3.second /= ftmp._plane3.second.norm();

      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane3;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.first = ctmp;
        // [point3, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.second = ctmp;
        computed[(*it)._point3][(*it)._point1] = computed[(*it)._point1][(*it)._point3] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane3;
        tVVR[tind].second.second.first = scount;
      }

      bsVVR.push_back(ftmp);
      ++scount;
    }

    //change the "outer BV" ids in the VVR limits to set the true values
    for(unsigned(i) = 0 ; i < ssVVR.size() ; ++i)
    {
      for(std::vector<SCD::VVRcone>::iterator it = ssVVR[i].begin() ; it != ssVVR[i].end() ; ++it)
      {
        it->first = smallSpheres.size() + bigSpheres.size() + it->first;
      }
    }
    for(unsigned(i) = 0 ; i < bsVVR.size() ; ++i)
    {
      bsVVR[i]._plane1.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first;
      bsVVR[i]._plane2.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first;
      bsVVR[i]._plane3.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first;
    }
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      it->second.second.first.first += smallSpheres.size();
      it->second.second.second.first += smallSpheres.size();
    }

    //sort the data to erase the useless toruses
    /*METHODE :
      - trouver les tores inutiles : ce sont ceux dont les deux grandes sphères de chaque côté ont le même centre
      - pour ça, parcourir la liste des VVR de tores et regarder les sphères de chaque côté.
      - noter les ID des tores à supprimer
      - pour chacun, noter les ID des BV avec lesquels il est en relation : petites sphères et grandes sphères
      - supprimer les tores inutiles
      - supprimer les VVR de ces tores
      - trouver les BV avec lesquels chacun était en relation
      - petites sphères : supprimer la VVR limite en question (point)
      - grande sphère : garder la VVR limite mais changer l'ID du outBV -> sphère de l'autre côté
      (en fait tout cette étape se fait uniquement au niveau des VVR des petites et grandes sphères, y'a rien qui change 
      au niveau des BV eux-mêmes)
     */
    std::map<int, bool> torusToErase, torusToEraseID;
    std::multimap<int, int> linkedsSphereIDs;
    std::multimap<int, std::pair<int, int> > linkedbSphereIDs;
    //find the useless toruses : those for which the two linked big spheres have the same center;
    //and store the data to make the needed changes after
    int currentTorusID;
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      int firstID = it->second.second.first.first - smallSpheres.size();
      int secondID = it->second.second.second.first - smallSpheres.size();
      if((bigSpheres[firstID]._center - bigSpheres[secondID]._center).normsquared() < epsilon)
      {
        currentTorusID = smallSpheres.size() + bigSpheres.size() + tinds[it->first];
        torusToEraseID.insert(std::pair<int, bool>(currentTorusID, true));
        torusToErase.insert(std::pair<int, bool>(it->first, true));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.first.first, currentTorusID));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.second.first, currentTorusID));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.first.first, std::pair<int, int>(currentTorusID, it->second.second.second.first)));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.second.first, std::pair<int, int>(currentTorusID, it->second.second.first.first)));
      }
    }

    //erase those toruses and their VVR limits
    /*std::map<int, bool>::iterator itm = torusToErase.end();
      --itm;
      for( ; itm != torusToErase.begin() ; --itm)
      {
      tVVR.erase(itm->first);

      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != itm->first)
      {
      --itt;
      }
      torus.erase(itt);
      }
      tVVR.erase(torusToErase.begin()->first);
      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != torusToErase.begin()->first)
      {
      --itt;
      }
      torus.erase(itt);*/

    //erase the VVR limits that separates them from the small spheres in the small spheres' VVR list
    for(std::multimap<int, int>::iterator it = linkedsSphereIDs.begin() ; it != linkedsSphereIDs.end() ; ++it)
    {
      std::vector<VVRcone>::iterator it2 = ssVVR[it->first].begin();
      while( (it2 != ssVVR[it->first].end()) && (it2->first != it->second) )
      {
        ++it2;
      }
      if(it2 != ssVVR[it->first].end())
        ssVVR[it->first].erase(it2);
    }

    //change the "outer BV" for the concerned big spheres' VVR limits : it's not the erased torus but the big 
    //sphere that has the same center
    for(std::multimap<int, std::pair<int, int> >::iterator it = linkedbSphereIDs.begin() ; it != linkedbSphereIDs.end() ; ++it)
    {
      int id = it->first - smallSpheres.size();
      if(bsVVR[id]._plane1.first == it->second.first)
        bsVVR[id]._plane1.first = it->second.second;
      else if(bsVVR[id]._plane2.first == it->second.first)
        bsVVR[id]._plane2.first = it->second.second;
      else if(bsVVR[id]._plane3.first == it->second.first)
        bsVVR[id]._plane3.first = it->second.second;
    }

    //change the "outer BV" IDs of the small and big spheres because some of them are too big
    for(std::map<int, bool>::reverse_iterator it = torusToEraseID.rbegin() ; it != torusToEraseID.rend() ; ++it)
    {
      for(std::vector< std::vector<VVRcone> >::iterator it2 = ssVVR.begin() ; it2 != ssVVR.end() ; ++it2)
      {
        for(std::vector<VVRcone>::iterator it3 = it2->begin() ; it3 != it2->end() ; ++it3)
        {
          if(it3->first > it->first)
            --(it3->first);
        }
      }
      for(std::vector<faceVVR>::iterator it2 = bsVVR.begin() ; it2 != bsVVR.end() ; ++it2)
      {
        if(it2->_plane1.first > it->first)
          --(it2->_plane1.first);
        if(it2->_plane2.first > it->first)
          --(it2->_plane2.first);
        if(it2->_plane3.first > it->first)
          --(it2->_plane3.first);
      }
    }



    //write data into a file
    os.open(filename.c_str());
    os.precision(16);
    os<<_r<<" "<<_R<<std::endl;

    os<<" "<< ind.size() << std::endl;
    i = 0;


    for(std::vector<int>::const_iterator it = ind.begin() ; it != ind.end() ; ++it)
    {
      os << smallSpheres[*it]._radius << " " << smallSpheres[*it]._center.x << " " << smallSpheres[*it]._center.y << " " << smallSpheres[*it]._center.z << std::endl;

      //VVR
      os << ssVVR[ind[i]].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[ind[i]].begin() ; it2 != ssVVR[ind[i]].end() ; ++it2)
      {
        os << it2->first-difference << " ";
        os << it2->second._cosangle << " ";
        os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }

      ++i;
    }
    i = 0;
    os << bigSpheres.size() << std::endl;
    for(std::vector<SCD::VVRsphere>::const_iterator it = bigSpheres.begin() ; it != bigSpheres.end() ; ++it)
    {
      os << it->_radius << " " << it->_center.x << " " << it->_center.y << " " << it->_center.z << std::endl;
      //face vertices
      os << _points[_spheres[i]._point1].x << " " << _points[_spheres[i]._point1].y << " " << _points[_spheres[i]._point1].z << " ";
      os << _points[_spheres[i]._point2].x << " " << _points[_spheres[i]._point2].y << " " << _points[_spheres[i]._point2].z << " ";
      os << _points[_spheres[i]._point3].x << " " << _points[_spheres[i]._point3].y << " " << _points[_spheres[i]._point3].z << std::endl;

      //VVR
      os << bsVVR[i]._plane1.first-difference << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << bsVVR[i]._plane2.first-difference << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << bsVVR[i]._plane3.first-difference << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;

      ++i;
    }
    os << torus.size() << std::endl;
    for(std::vector< std::pair<int, VVRtorus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
    {
      if(torusToErase.find(it->first) != torusToErase.end())
        os << false << std::endl;
      else
        os << true << std::endl;
      os << it->second._extRadius << " ";
      os << _R << " ";
      os << it->second._center.x << " " << it->second._center.y << " " << it->second._center.z << " ";
      os << it->second._normal.x << " " << it->second._normal.y << " " << it->second._normal.z << std::endl;

      //VVR
      os << invind[tVVR[it->first].first.first.first] << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << invind[tVVR[it->first].first.second.first] << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << tVVR[it->first].second.first.first-difference << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << tVVR[it->first].second.second.first-difference << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
    }

    /*for(int i = 0 ; i < smallSpheres.size() ; ++i)
      {
      os << ssVVR[i].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[i].begin() ; it2 != ssVVR[i].end() ; ++it2)
      {
      os << (smallSpheres.size() + bigSpheres.size() + it2->first) << " ";
      os << it2->second._cosangle << " ";
      os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }
      }
      for(int i = 0; i < bigSpheres.size() ; ++i)
      {
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first) << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first) << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first) << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;
      }
      for(std::vector< std::pair<int, Torus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
      {
      os << tVVR[it->first].first.first.first << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << tVVR[it->first].first.second.first << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.first.first) << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.second.first) << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
      }*/
    os.close();
  }

  //#ifdef FALSE
  void SmoothHullGeneratorVVR::computeVVR(const std::string& filename)
  {
    std::ofstream os;
    std::vector<VVRsphere> smallSpheres;
    std::vector<VVRsphere> bigSpheres;
    std::vector< std::pair<int, VVRtorus> > torus;
    std::map<int, int> tinds;
    std::vector<std::vector<VVRcone> > ssVVR;
    std::vector<std::vector<bool> > computed(_points.size());
    std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > > tVVR;
    std::vector<SCD::faceVVR> bsVVR;
    VVRsphere stmp;
    VVRcone ctmp;
    double l;
    VVRtorus ttmp;
    int scount = 0;
    int tind, i;
    std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > tpair;
    SCD::faceVVR ftmp;
    double epsilon = 1e-10;

    //computing hull
    _index.clear();
    _spheres.clear();
    cover();

    for(i = 0 ; i < _points.size() ; ++i)
    {
      computed[i].resize(_points.size());
      for(int j = 0 ; j < _points.size() ; ++j)
        computed[i][j] = false;
    }

    for(std::vector<vector3d>::const_iterator it = _points.begin() ;
        it != _points.end() ;
        ++it)
    {
      stmp._center = *it;
      stmp._radius = _r;
      smallSpheres.push_back(stmp);
    }
    ssVVR.resize(_points.size());

    for(std::vector<VVRSFace>::const_iterator it = _spheres.begin() ;
        it != _spheres.end() ;
        ++it)
    {
      //register the new big sphere
      stmp._center = (*it)._center;
      stmp._radius = _R;
      bigSpheres.push_back(stmp);

      tind = getKey((*it)._point1, (*it)._point2);
      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point1] + _points[(*it)._point2]) / 2;
        l = (_points[(*it)._point1] - _points[(*it)._point2]).norm();
        ttmp._normal = (_points[(*it)._point1] - _points[(*it)._point2]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane1.first = tinds[tind];
      ftmp._plane1.second = (_points[(*it)._point1] - stmp._center)^(_points[(*it)._point2] - stmp._center);
      ftmp._plane1.second /= ftmp._plane1.second.norm();

      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the VVR plane for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane1;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.first = ctmp;
        // [point2, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.second = ctmp;
        computed[(*it)._point1][(*it)._point2] = computed[(*it)._point2][(*it)._point1] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR plane for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane1;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point2, (*it)._point3);
      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point2] + _points[(*it)._point3]) / 2;
        l = (_points[(*it)._point2] - _points[(*it)._point3]).norm();
        ttmp._normal = (_points[(*it)._point2] - _points[(*it)._point3]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane2.first = tinds[tind];
      ftmp._plane2.second = (_points[(*it)._point2] - stmp._center)^(_points[(*it)._point3] - stmp._center);
      ftmp._plane2.second /= ftmp._plane2.second.norm();

      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane2;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point2, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.first = ctmp;
        // [point3, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.second = ctmp;
        computed[(*it)._point2][(*it)._point3] = computed[(*it)._point3][(*it)._point2] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane2;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point3, (*it)._point1);
      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point3] + _points[(*it)._point1]) / 2;
        l = (_points[(*it)._point3] - _points[(*it)._point1]).norm();
        ttmp._normal = (_points[(*it)._point3] - _points[(*it)._point1]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane3.first = tinds[tind];
      ftmp._plane3.second = (_points[(*it)._point3] - stmp._center)^(_points[(*it)._point1] - stmp._center);
      ftmp._plane3.second /= ftmp._plane3.second.norm();

      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane3;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.first = ctmp;
        // [point3, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.second = ctmp;
        computed[(*it)._point3][(*it)._point1] = computed[(*it)._point1][(*it)._point3] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane3;
        tVVR[tind].second.second.first = scount;
      }

      bsVVR.push_back(ftmp);
      ++scount;
    }

    //change the "outer BV" ids in the VVR limits to set the true values
    for(i = 0 ; i < ssVVR.size() ; ++i)
    {
      for(std::vector<SCD::VVRcone>::iterator it = ssVVR[i].begin() ; it != ssVVR[i].end() ; ++it)
      {
        it->first = smallSpheres.size() + bigSpheres.size() + it->first;
      }
    }
    for(i = 0 ; i < bsVVR.size() ; ++i)
    {
      bsVVR[i]._plane1.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first;
      bsVVR[i]._plane2.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first;
      bsVVR[i]._plane3.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first;
    }
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      it->second.second.first.first += smallSpheres.size();
      it->second.second.second.first += smallSpheres.size();
    }

    //sort the data to erase the useless toruses
    /*METHODE :
      - trouver les tores inutiles : ce sont ceux dont les deux grandes sphères de chaque côté ont le même centre
      - pour ça, parcourir la liste des VVR de tores et regarder les sphères de chaque côté.
      - noter les ID des tores à supprimer
      - pour chacun, noter les ID des BV avec lesquels il est en relation : petites sphères et grandes sphères
      - supprimer les tores inutiles
      - supprimer les VVR de ces tores
      - trouver les BV avec lesquels chacun était en relation
      - petites sphères : supprimer la VVR limite en question (point)
      - grande sphère : garder la VVR limite mais changer l'ID du outBV -> sphère de l'autre côté
      (en fait tout cette étape se fait uniquement au niveau des VVR des petites et grandes sphères, y'a rien qui change
      au niveau des BV eux-mêmes)
     */
    std::map<int, bool> torusToErase, torusToEraseID;
    std::multimap<int, int> linkedsSphereIDs;
    std::multimap<int, std::pair<int, int> > linkedbSphereIDs;
    //find the useless toruses : those for which the two linked big spheres have the same center;
    //and store the data to make the needed changes after
    int currentTorusID;
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      int firstID = it->second.second.first.first - smallSpheres.size();
      int secondID = it->second.second.second.first - smallSpheres.size();
      if((bigSpheres[firstID]._center - bigSpheres[secondID]._center).normsquared() < epsilon)
      {
        currentTorusID = smallSpheres.size() + bigSpheres.size() + tinds[it->first];
        torusToEraseID.insert(std::pair<int, bool>(currentTorusID, true));
        torusToErase.insert(std::pair<int, bool>(it->first, true));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.first.first, currentTorusID));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.second.first, currentTorusID));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.first.first, std::pair<int, int>(currentTorusID, it->second.second.second.first)));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.second.first, std::pair<int, int>(currentTorusID, it->second.second.first.first)));
      }
    }

    //erase those toruses and their VVR limits
    /*std::map<int, bool>::iterator itm = torusToErase.end();
      --itm;
      for( ; itm != torusToErase.begin() ; --itm)
      {
      tVVR.erase(itm->first);

      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != itm->first)
      {
      --itt;
      }
      torus.erase(itt);
      }
      tVVR.erase(torusToErase.begin()->first);
      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != torusToErase.begin()->first)
      {
      --itt;
      }
      torus.erase(itt);*/

    //erase the VVR limits that separates them from the small spheres in the small spheres' VVR list
    for(std::multimap<int, int>::iterator it = linkedsSphereIDs.begin() ; it != linkedsSphereIDs.end() ; ++it)
    {
      std::vector<VVRcone>::iterator it2 = ssVVR[it->first].begin();
      while( (it2 != ssVVR[it->first].end()) && (it2->first != it->second) )
      {
        ++it2;
      }
      if(it2 != ssVVR[it->first].end())
        ssVVR[it->first].erase(it2);
    }

    //change the "outer BV" for the concerned big spheres' VVR limits : it's not the erased torus but the big
    //sphere that has the same center
    for(std::multimap<int, std::pair<int, int> >::iterator it = linkedbSphereIDs.begin() ; it != linkedbSphereIDs.end() ; ++it)
    {
      int id = it->first - smallSpheres.size();
      if(bsVVR[id]._plane1.first == it->second.first)
        bsVVR[id]._plane1.first = it->second.second;
      else if(bsVVR[id]._plane2.first == it->second.first)
        bsVVR[id]._plane2.first = it->second.second;
      else if(bsVVR[id]._plane3.first == it->second.first)
        bsVVR[id]._plane3.first = it->second.second;
    }

    //change the "outer BV" IDs of the small and big spheres because some of them are too big
    for(std::map<int, bool>::reverse_iterator it = torusToEraseID.rbegin() ; it != torusToEraseID.rend() ; ++it)
    {
      for(std::vector< std::vector<VVRcone> >::iterator it2 = ssVVR.begin() ; it2 != ssVVR.end() ; ++it2)
      {
        for(std::vector<VVRcone>::iterator it3 = it2->begin() ; it3 != it2->end() ; ++it3)
        {
          if(it3->first > it->first)
            --(it3->first);
        }
      }
      for(std::vector<faceVVR>::iterator it2 = bsVVR.begin() ; it2 != bsVVR.end() ; ++it2)
      {
        if(it2->_plane1.first > it->first)
          --(it2->_plane1.first);
        if(it2->_plane2.first > it->first)
          --(it2->_plane2.first);
        if(it2->_plane3.first > it->first)
          --(it2->_plane3.first);
      }
    }



    //write data into a file
    os.open(filename.c_str());
    os.precision(16);
    i = 0;


    os<<_r<<" "<<_R<<std::endl;

    os << smallSpheres.size() << std::endl;
    for(std::vector<SCD::VVRsphere>::const_iterator it = smallSpheres.begin() ; it != smallSpheres.end() ; ++it)
    {
      os << it->_radius << " " << it->_center.x << " " << it->_center.y << " " << it->_center.z << std::endl;

      //VVR
      os << ssVVR[i].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[i].begin() ; it2 != ssVVR[i].end() ; ++it2)
      {
        os << it2->first << " ";
        os << it2->second._cosangle << " ";
        os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }

      ++i;
    }
    i = 0;
    os << bigSpheres.size() << std::endl;
    for(std::vector<SCD::VVRsphere>::const_iterator it = bigSpheres.begin() ; it != bigSpheres.end() ; ++it)
    {
      os << it->_radius << " " << it->_center.x << " " << it->_center.y << " " << it->_center.z << std::endl;
      //face vertices
      os << _points[_spheres[i]._point1].x << " " << _points[_spheres[i]._point1].y << " " << _points[_spheres[i]._point1].z << " ";
      os << _points[_spheres[i]._point2].x << " " << _points[_spheres[i]._point2].y << " " << _points[_spheres[i]._point2].z << " ";
      os << _points[_spheres[i]._point3].x << " " << _points[_spheres[i]._point3].y << " " << _points[_spheres[i]._point3].z << std::endl;

      //VVR
      os << bsVVR[i]._plane1.first << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << bsVVR[i]._plane2.first << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << bsVVR[i]._plane3.first << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;

      ++i;
    }
    os << torus.size() << std::endl;
    for(std::vector< std::pair<int, VVRtorus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
    {
      if(torusToErase.find(it->first) != torusToErase.end())
        os << false << std::endl;
      else
        os << true << std::endl;
      os << it->second._extRadius << " ";
      os << _R << " ";
      os << it->second._center.x << " " << it->second._center.y << " " << it->second._center.z << " ";
      os << it->second._normal.x << " " << it->second._normal.y << " " << it->second._normal.z << std::endl;

      //VVR
      os << tVVR[it->first].first.first.first << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << tVVR[it->first].first.second.first << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << tVVR[it->first].second.first.first << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << tVVR[it->first].second.second.first << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
    }

    /*for(int i = 0 ; i < smallSpheres.size() ; ++i)
      {
      os << ssVVR[i].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[i].begin() ; it2 != ssVVR[i].end() ; ++it2)
      {
      os << (smallSpheres.size() + bigSpheres.size() + it2->first) << " ";
      os << it2->second._cosangle << " ";
      os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }
      }
      for(int i = 0; i < bigSpheres.size() ; ++i)
      {
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first) << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first) << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first) << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;
      }
      for(std::vector< std::pair<int, Torus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
      {
      os << tVVR[it->first].first.first.first << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << tVVR[it->first].first.second.first << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.first.first) << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.second.first) << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
      }*/
    os.close();
  }
  //#endif
  void SmoothHullGeneratorVVR::computeVVR_WithPolyhedron(const std::string& filename)
  {
    std::ofstream os;
    std::vector<VVRsphere> smallSpheres;
    std::vector<VVRsphere> bigSpheres;
    std::vector< std::pair<int, VVRtorus> > torus;
    std::map<int, int> tinds;
    std::vector<std::vector<VVRcone> > ssVVR;
    std::vector<std::vector<bool> > computed(_points.size());
    std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > > tVVR;
    std::vector<SCD::faceVVR> bsVVR;
    VVRsphere stmp;
    VVRcone ctmp;
    double l;
    VVRtorus ttmp;
    int scount = 0;
    int tind, i;
    std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > tpair;
    SCD::faceVVR ftmp;
    double epsilon = 1e-10;

    //computing hull
    _index.clear();
    _spheres.clear();
    cover();

    std::vector<int> ind(_index.size());
    std::vector<int> invind(_points.size(),-1);

    int difference=_points.size()-_index.size();

    i=0;
    for (std::set<int,std::less<int> >::iterator it=_index.begin(); it!=_index.end() ; ++it,++i)
    {
      ind[i]=*it;
      invind[*it]=i;
    }

    for(i = 0 ; i < _points.size() ; ++i)
    {
      computed[i].resize(_points.size());
      for(int j = 0 ; j < _points.size() ; ++j)
        computed[i][j] = false;
    }

    for(std::vector<vector3d>::const_iterator it = _points.begin() ;
        it != _points.end() ;
        ++it)
    {
      stmp._center = *it;
      stmp._radius = _r;
      smallSpheres.push_back(stmp);
    }
    ssVVR.resize(_points.size());

    for(std::vector<VVRSFace>::const_iterator it = _spheres.begin() ;
        it != _spheres.end() ;
        ++it)
    {
      //register the new big sphere
      stmp._center = (*it)._center;
      stmp._radius = _R;
      bigSpheres.push_back(stmp);

      tind = getKey((*it)._point1, (*it)._point2);
      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point1] + _points[(*it)._point2]) / 2;
        l = (_points[(*it)._point1] - _points[(*it)._point2]).norm();
        ttmp._normal = (_points[(*it)._point1] - _points[(*it)._point2]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane1.first = tinds[tind];
      ftmp._plane1.second = (_points[(*it)._point1] - stmp._center)^(_points[(*it)._point2] - stmp._center);
      ftmp._plane1.second /= ftmp._plane1.second.norm();

      if(!computed[(*it)._point1][(*it)._point2])
      {
        //register the VVR plane for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane1;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.first = ctmp;
        // [point2, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.second = ctmp;
        computed[(*it)._point1][(*it)._point2] = computed[(*it)._point2][(*it)._point1] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR plane for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane1;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point2, (*it)._point3);
      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point2] + _points[(*it)._point3]) / 2;
        l = (_points[(*it)._point2] - _points[(*it)._point3]).norm();
        ttmp._normal = (_points[(*it)._point2] - _points[(*it)._point3]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane2.first = tinds[tind];
      ftmp._plane2.second = (_points[(*it)._point2] - stmp._center)^(_points[(*it)._point3] - stmp._center);
      ftmp._plane2.second /= ftmp._plane2.second.norm();

      if(!computed[(*it)._point2][(*it)._point3])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane2;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point2, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point2].push_back(ctmp);
        ctmp.first = (*it)._point2;
        tpair.first.first = ctmp;
        // [point3, point2]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.second = ctmp;
        computed[(*it)._point2][(*it)._point3] = computed[(*it)._point3][(*it)._point2] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane2;
        tVVR[tind].second.second.first = scount;
      }

      tind = getKey((*it)._point3, (*it)._point1);
      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the torus
        ttmp._center = (_points[(*it)._point3] + _points[(*it)._point1]) / 2;
        l = (_points[(*it)._point3] - _points[(*it)._point1]).norm();
        ttmp._normal = (_points[(*it)._point3] - _points[(*it)._point1]) / l;
        ttmp._extRadius = (ttmp._center - stmp._center).norm();
        torus.push_back(std::pair<int, VVRtorus>(tind, ttmp));
        tinds.insert(std::pair<int, int>(tind, torus.size() - 1));
      }

      ftmp._plane3.first = tinds[tind];
      ftmp._plane3.second = (_points[(*it)._point3] - stmp._center)^(_points[(*it)._point1] - stmp._center);
      ftmp._plane3.second /= ftmp._plane3.second.norm();

      if(!computed[(*it)._point3][(*it)._point1])
      {
        //register the VVR cone for the big sphere and torus in the torus
        tpair.second.first = ftmp._plane3;
        tpair.second.first.first = scount;

        //register the VVR cone for the small sphere and torus
        // [point1, point3]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = ttmp._normal;
        ctmp.second._cosangle = l / ( 2 * (_R - _r));
        ssVVR[(*it)._point3].push_back(ctmp);
        ctmp.first = (*it)._point3;
        tpair.first.first = ctmp;
        // [point3, point1]
        ctmp.first = torus.size() - 1;
        ctmp.second._axis = -ctmp.second._axis;
        ssVVR[(*it)._point1].push_back(ctmp);
        ctmp.first = (*it)._point1;
        tpair.first.second = ctmp;
        computed[(*it)._point3][(*it)._point1] = computed[(*it)._point1][(*it)._point3] = true;

        tVVR.insert(std::pair<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >(tind, tpair));
      }
      else
      {
        //register the VVR cone for the big sphere and torus
        tVVR[tind].second.second = ftmp._plane3;
        tVVR[tind].second.second.first = scount;
      }

      bsVVR.push_back(ftmp);
      ++scount;
    }

    //change the "outer BV" ids in the VVR limits to set the true values
    for(i = 0 ; i < ssVVR.size() ; ++i)
    {
      for(std::vector<SCD::VVRcone>::iterator it = ssVVR[i].begin() ; it != ssVVR[i].end() ; ++it)
      {
        it->first = smallSpheres.size() + bigSpheres.size() + it->first;
      }
    }
    for(i = 0 ; i < bsVVR.size() ; ++i)
    {
      bsVVR[i]._plane1.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first;
      bsVVR[i]._plane2.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first;
      bsVVR[i]._plane3.first = smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first;
    }
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      it->second.second.first.first += smallSpheres.size();
      it->second.second.second.first += smallSpheres.size();
    }

    //sort the data to erase the useless toruses
    /*METHODE :
      - trouver les tores inutiles : ce sont ceux dont les deux grandes sphères de chaque côté ont le même centre
      - pour ça, parcourir la liste des VVR de tores et regarder les sphères de chaque côté.
      - noter les ID des tores à supprimer
      - pour chacun, noter les ID des BV avec lesquels il est en relation : petites sphères et grandes sphères
      - supprimer les tores inutiles
      - supprimer les VVR de ces tores
      - trouver les BV avec lesquels chacun était en relation
      - petites sphères : supprimer la VVR limite en question (point)
      - grande sphère : garder la VVR limite mais changer l'ID du outBV -> sphère de l'autre côté
      (en fait tout cette étape se fait uniquement au niveau des VVR des petites et grandes sphères, y'a rien qui change
      au niveau des BV eux-mêmes)
     */
    std::map<int, bool> torusToErase, torusToEraseID;
    std::multimap<int, int> linkedsSphereIDs;
    std::multimap<int, std::pair<int, int> > linkedbSphereIDs;
    //find the useless toruses : those for which the two linked big spheres have the same center;
    //and store the data to make the needed changes after
    int currentTorusID;
    for(std::map<int, std::pair<std::pair<VVRcone, VVRcone>, std::pair<VVRplane, VVRplane> > >::iterator it = tVVR.begin() ;
        it != tVVR.end() ; ++it)
    {
      int firstID = it->second.second.first.first - smallSpheres.size();
      int secondID = it->second.second.second.first - smallSpheres.size();
      if((bigSpheres[firstID]._center - bigSpheres[secondID]._center).normsquared() < epsilon)
      {
        currentTorusID = smallSpheres.size() + bigSpheres.size() + tinds[it->first];
        torusToEraseID.insert(std::pair<int, bool>(currentTorusID, true));
        torusToErase.insert(std::pair<int, bool>(it->first, true));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.first.first, currentTorusID));
        linkedsSphereIDs.insert(std::pair<int, int>(it->second.first.second.first, currentTorusID));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.first.first, std::pair<int, int>(currentTorusID, it->second.second.second.first)));
        linkedbSphereIDs.insert(std::pair<int, std::pair<int, int> >(it->second.second.second.first, std::pair<int, int>(currentTorusID, it->second.second.first.first)));
      }
    }

    //erase those toruses and their VVR limits
    /*std::map<int, bool>::iterator itm = torusToErase.end();
      --itm;
      for( ; itm != torusToErase.begin() ; --itm)
      {
      tVVR.erase(itm->first);

      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != itm->first)
      {
      --itt;
      }
      torus.erase(itt);
      }
      tVVR.erase(torusToErase.begin()->first);
      std::vector< std::pair<int, Torus> >::iterator itt = torus.end();
      --itt;
      while(itt->first != torusToErase.begin()->first)
      {
      --itt;
      }
      torus.erase(itt);*/

    //erase the VVR limits that separates them from the small spheres in the small spheres' VVR list
    for(std::multimap<int, int>::iterator it = linkedsSphereIDs.begin() ; it != linkedsSphereIDs.end() ; ++it)
    {
      std::vector<VVRcone>::iterator it2 = ssVVR[it->first].begin();
      while( (it2 != ssVVR[it->first].end()) && (it2->first != it->second) )
      {
        ++it2;
      }
      if(it2 != ssVVR[it->first].end())
        ssVVR[it->first].erase(it2);
    }

    //change the "outer BV" for the concerned big spheres' VVR limits : it's not the erased torus but the big
    //sphere that has the same center
    for(std::multimap<int, std::pair<int, int> >::iterator it = linkedbSphereIDs.begin() ; it != linkedbSphereIDs.end() ; ++it)
    {
      int id = it->first - smallSpheres.size();
      if(bsVVR[id]._plane1.first == it->second.first)
        bsVVR[id]._plane1.first = it->second.second;
      else if(bsVVR[id]._plane2.first == it->second.first)
        bsVVR[id]._plane2.first = it->second.second;
      else if(bsVVR[id]._plane3.first == it->second.first)
        bsVVR[id]._plane3.first = it->second.second;
    }

    //change the "outer BV" IDs of the small and big spheres because some of them are too big
    for(std::map<int, bool>::reverse_iterator it = torusToEraseID.rbegin() ; it != torusToEraseID.rend() ; ++it)
    {
      for(std::vector< std::vector<VVRcone> >::iterator it2 = ssVVR.begin() ; it2 != ssVVR.end() ; ++it2)
      {
        for(std::vector<VVRcone>::iterator it3 = it2->begin() ; it3 != it2->end() ; ++it3)
        {
          if(it3->first > it->first)
            --(it3->first);
        }
      }
      for(std::vector<faceVVR>::iterator it2 = bsVVR.begin() ; it2 != bsVVR.end() ; ++it2)
      {
        if(it2->_plane1.first > it->first)
          --(it2->_plane1.first);
        if(it2->_plane2.first > it->first)
          --(it2->_plane2.first);
        if(it2->_plane3.first > it->first)
          --(it2->_plane3.first);
      }
    }



    //write data into a file
    os.open(filename.c_str());
    os.precision(16);
    os<<_r<<" "<<_R<<std::endl;

    os<<" "<< ind.size() << std::endl;
    i = 0;


    for(std::vector<int>::const_iterator it = ind.begin() ; it != ind.end() ; ++it)
    {
      os << smallSpheres[*it]._radius << " " << smallSpheres[*it]._center.x << " " << smallSpheres[*it]._center.y << " " << smallSpheres[*it]._center.z << std::endl;

      //VVR
      os << ssVVR[ind[i]].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[ind[i]].begin() ; it2 != ssVVR[ind[i]].end() ; ++it2)
      {
        os << it2->first-difference << " ";
        os << it2->second._cosangle << " ";
        os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }

      ++i;
    }
    i = 0;
    os << bigSpheres.size() << std::endl;
    for(std::vector<SCD::VVRsphere>::const_iterator it = bigSpheres.begin() ; it != bigSpheres.end() ; ++it)
    {
      os << it->_radius << " " << it->_center.x << " " << it->_center.y << " " << it->_center.z << std::endl;
      //face vertices
      os << _points[_spheres[i]._point1].x << " " << _points[_spheres[i]._point1].y << " " << _points[_spheres[i]._point1].z << " ";
      os << _points[_spheres[i]._point2].x << " " << _points[_spheres[i]._point2].y << " " << _points[_spheres[i]._point2].z << " ";
      os << _points[_spheres[i]._point3].x << " " << _points[_spheres[i]._point3].y << " " << _points[_spheres[i]._point3].z << std::endl;

      //VVR
      os << bsVVR[i]._plane1.first-difference << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << bsVVR[i]._plane2.first-difference << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << bsVVR[i]._plane3.first-difference << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;

      ++i;
    }
    os << torus.size() << std::endl;
    for(std::vector< std::pair<int, VVRtorus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
    {
      if(torusToErase.find(it->first) != torusToErase.end())
        os << false << std::endl;
      else
        os << true << std::endl;
      os << it->second._extRadius << " ";
      os << _R << " ";
      os << it->second._center.x << " " << it->second._center.y << " " << it->second._center.z << " ";
      os << it->second._normal.x << " " << it->second._normal.y << " " << it->second._normal.z << std::endl;

      //VVR
      os << invind[tVVR[it->first].first.first.first] << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << invind[tVVR[it->first].first.second.first] << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << tVVR[it->first].second.first.first-difference << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << tVVR[it->first].second.second.first-difference << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
    }

    /*for(int i = 0 ; i < smallSpheres.size() ; ++i)
      {
      os << ssVVR[i].size() << std::endl;
      for(std::vector<SCD::VVRcone>::const_iterator it2 = ssVVR[i].begin() ; it2 != ssVVR[i].end() ; ++it2)
      {
      os << (smallSpheres.size() + bigSpheres.size() + it2->first) << " ";
      os << it2->second._cosangle << " ";
      os << it2->second._axis.x << " " << it2->second._axis.y << " " << it2->second._axis.z << std::endl;
      }
      }
      for(int i = 0; i < bigSpheres.size() ; ++i)
      {
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane1.first) << " ";
      os << bsVVR[i]._plane1.second.x << " " << bsVVR[i]._plane1.second.y << " " << bsVVR[i]._plane1.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane2.first) << " ";
      os << bsVVR[i]._plane2.second.x << " " << bsVVR[i]._plane2.second.y << " " << bsVVR[i]._plane2.second.z << std::endl;
      os << (smallSpheres.size() + bigSpheres.size() + bsVVR[i]._plane3.first) << " ";
      os << bsVVR[i]._plane3.second.x << " " << bsVVR[i]._plane3.second.y << " " << bsVVR[i]._plane3.second.z << std::endl;
      }
      for(std::vector< std::pair<int, Torus> >::const_iterator it = torus.begin() ; it != torus.end() ; ++it)
      {
      os << tVVR[it->first].first.first.first << " ";
      os << tVVR[it->first].first.first.second._cosangle << " ";
      os << tVVR[it->first].first.first.second._axis.x << " ";
      os << tVVR[it->first].first.first.second._axis.y << " ";
      os << tVVR[it->first].first.first.second._axis.z << std::endl;
      os << tVVR[it->first].first.second.first << " ";
      os << tVVR[it->first].first.second.second._cosangle << " ";
      os << tVVR[it->first].first.second.second._axis.x << " ";
      os << tVVR[it->first].first.second.second._axis.y << " ";
      os << tVVR[it->first].first.second.second._axis.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.first.first) << " ";
      os << tVVR[it->first].second.first.second.x << " ";
      os << tVVR[it->first].second.first.second.y << " ";
      os << tVVR[it->first].second.first.second.z << std::endl;
      os << (smallSpheres.size() + tVVR[it->first].second.second.first) << " ";
      os << tVVR[it->first].second.second.second.x << " ";
      os << tVVR[it->first].second.second.second.y << " ";
      os << tVVR[it->first].second.second.second.z << std::endl;
      }*/
    os.close();


    os.open((filename+std::string(".inp")).c_str());
    os.precision(16);
    i = 0;


    os<<"3 ";

    os <<  ind.size() << std::endl;
    for(std::vector<int>::const_iterator it = ind.begin() ; it != ind.end() ; ++it)
    {
      os  << smallSpheres[*it]._center.x << " " << smallSpheres[*it]._center.y << " " << smallSpheres[*it]._center.z << std::endl;
    }
    os.close();


    //spawnl(P_WAIT,"qconvex.exe","qconvex.exe",(std::string(" TI ")+filename+std::string(".inp")).c_str(),(std::string(" TO ")+filename+std::string(".otp")).c_str()," Qt"," o"," f",NULL);//appel à qconvex.exe


  }


}
