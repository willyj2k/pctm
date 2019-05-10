//
// Created by Rene Lee on 2019-04-30.
//

#ifndef PCTM_POINT_H
#define PCTM_POINT_H

#include "CGL/CGL.h"

using namespace std;
using namespace CGL;

class Point {
public:
  Vector3D pos;
  Vector3D normal;

  /**
   * Constructor.
   */
  Point(Vector3D p, Vector3D n) : pos( p ), normal( n ) { }
  Point(Vector3D p) : pos( p ), normal( Vector3D(0, 0, 0) ) { }
};


#endif //PCTM_POINT_H
