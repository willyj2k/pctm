//
// Created by Rene Lee on 2019-04-23.
//

#ifndef BALLPIVOT_H
#define BALLPIVOT_H

#include "CGL/CGL.h"
#include "point.h"
#include <unordered_map>
#include <unordered_set>

using namespace CGL;

class BallPivot {
  public:
    void init(const std::vector<Point> &points, double radius, Vector3D bound_min, Vector3D bound_max);
    std::vector<Point *> find_seed_triangle();
    static double dist(const Point &p);

  private:
    // cell info for hashing
    double cell_width;

    struct CellIndex {
      int x_ind;
      int y_ind;
      int z_ind;

      // constructors
      CellIndex(int x, int y, int z) : x_ind( x ), y_ind( y ), z_ind( z ) { }

      // default constructor; be careful with this
      CellIndex() : x_ind( 0 ), y_ind( 0 ), z_ind( 0 ) { }

      // equality override
      bool operator==(const CellIndex &other) const {
        return (x_ind == other.x_ind
                && y_ind == other.y_ind
                && z_ind == other.z_ind);
      }

      bool operator!=(const CellIndex &other) const {
        return (x_ind != other.x_ind
                || y_ind == other.y_ind
                || z_ind == other.z_ind);
      }
    };

    // spatial map of points in the cloud
    //
    // use this as the master data structure for storing points;
    // we don't want to use pointers to Points here anymore and
    // instead want to just store them here raw
    std::unordered_map<int, std::vector<Point> *> spatial_map;

    // radius of the ball to be pivoted
    double radius;

    // vectors to define the bounding box of the point cloud for hashing
    // TODO not sure if we really need bound_max for now
    Vector3D bound_min;
    Vector3D bound_max;
  
    // set to track cells that have already been seeded and pivoted
    std::unordered_set<int> processed_cells;
    CellIndex seed_cell;
    CellIndex max_cell;

    // primes for hashing
    int large_prime = 2038074743;
    int small_prime = 113;

    void create_spatial_grid(const vector<Point> &points);
    vector<Point *> neighborhood(double r, const Point &p);
    Vector3D circumcenter(const Point &a, const Point &b, const Point &c);
    Vector3D ball_center(const Point &a, const Point &b, const Point &c, const Vector3D &normal);
    Vector3D naive_plane_normal(const Point &a, const Point &b, const Point &c);
    Vector3D correct_plane_normal(const Point &a, const Point &b, const Point &c);
    int hash_position(const Point &p);
    int hash_cell(const CellIndex &c);
    CellIndex get_cell(const Point &p);
    Point* get_seed_candidate(const CellIndex &c);
    void increment_seed_cell();
    void calculate_normals();
};


#endif //BALLPIVOT_H
