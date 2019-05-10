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
    void init(std::vector<Point> points, double radius, Vector3D bound_min);
    std::vector<Point> find_seed_triangle();
    static double dist(const Point &p);

  private:
    // std::vector of used points
    std::vector<Point> used;

    // std::vector of unused points
    std::vector<Point> unused;

    // spatial map
    std::unordered_map<double, std::vector<Point *> *> map;

    // radius of the ball to be pivoted
    double radius;

    // vector to define the minimum corner of the
    // bounding box of the point cloud for hashing
    Vector3D bound_min;
  
    // primes for hashing
    int large_prime = 29996224275833;
    int small_prime = 113;

    // cell width for hashing
    double cell_width = 2 * radius;

    void create_spatial_grid();
    Vector3D circumcenter(const Point &a, const Point &b, const Point &c);
    Vector3D ball_center(const Point &a, const Point &b, const Point &c, const Vector3D &normal);
    Vector3D naive_plane_normal(const Point &a, const Point &b, const Point &c);
    Vector3D correct_plane_normal(const Point &a, const Point &b, const Point &c);
    int hash_position(const Point &p);
    void calculate_normals();
};


#endif //BALLPIVOT_H
