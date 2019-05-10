//
// Created by Rene Lee on 2019-04-23.
//

#ifndef BALLPIVOT_H
#define BALLPIVOT_H

#include "CGL/CGL.h"
#include "point.h"
#include <unordered_map>

using namespace std;
using namespace CGL;

class BallPivot {
  public:
    void init(std::vector<Point> points, float radius);
    std::vector<Point> find_seed_triangle();

  private:
    // std::vector of used points
    std::vector<Point> used;

    // std::vector of unused points
    std::vector<Point> unused;

    // spatial map
    std::unordered_map<float, vector<Point *> *> map;

    // float radius
    float radius;

    double width;
    double height;

    Point *point;

    void create_spatial_grid ();
//    std::vector<Point *> find_seed_triangle();
    Vector3D circumcenter(const Point &a, const Point &b, const Point &c);
    Vector3D rho_center(double rho, const Point &a, const Point &b, const Point &c);
    Vector3D naive_plane_normal(const Point &a, const Point &b, const Point &c);
    Vector3D correct_plane_normal(const Point &a, const Point &b, const Point &c);
    float hash_position(const Point &p);
    float distance(const Point &a, const Point &b);
    bool compare(Point *a, Point *b);
};


#endif //BALLPIVOT_H
