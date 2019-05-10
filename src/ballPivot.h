//
// Created by Rene Lee on 2019-04-23.
//

#ifndef BALLPIVOT_H
#define BALLPIVOT_H

#include "CGL/CGL.h"
#include "point.h"
#include <unordered_map>

using namespace CGL;

class BallPivot {
  public:
    void init(std::vector<Point> points, double radius);
    std::vector<Point> find_seed_triangle();
    std::vector<Point> all_points;
    static double dist(const Point &p);

  private:
    // std::vector of used points
    std::vector<Point> used;

    // std::vector of unused points
    std::vector<Point> unused;

    // spatial map
    std::unordered_map<double, std::vector<Point *> *> map;

    double radius;
    double width;
    double height;

    void create_spatial_grid ();
    Vector3D circumcenter(const Point &a, const Point &b, const Point &c);
    Vector3D ball_center(const Point &a, const Point &b, const Point &c);
    Vector3D naive_plane_normal(const Point &a, const Point &b, const Point &c);
    Vector3D correct_plane_normal(const Point &a, const Point &b, const Point &c);
    double hash_position(const Point &p);
    void calculate_normals();
};


#endif //BALLPIVOT_H
