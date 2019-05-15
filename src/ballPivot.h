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
    // info for storing vertex information relating to ball pivoting
    enum class VertexSpecifier {i, j, k};

    struct PivotTriangle {
      Point *sigma_i;
      Point *sigma_j;
      Point *sigma_o;
      Point *center;
      bool empty;
      bool isBoundary;
      bool isFrozen;

      // constructors
      PivotTriangle(Point *i, Point *j, Point *o, Point *center) :
        sigma_i( i ),
        sigma_j( j ),
        sigma_o( o ),
        center( center ),
        empty( false ),
        isBoundary( false ),
        isFrozen( false )
        { }

      PivotTriangle() :
        sigma_i( NULL ),
        sigma_j( NULL ),
        sigma_o( NULL ),
        center( NULL ),
        empty( true ),
        isBoundary( false ),
        isFrozen( false )
        { }

      PivotTriangle(const PivotTriangle &other) :
        sigma_i( other.sigma_i ),
        sigma_j( other.sigma_j ),
        sigma_o( other.sigma_o ),
        center( other.center ),
        empty( other.empty ),
        isBoundary( other.isBoundary ),
        isFrozen( other.isFrozen )
        { }
    };

    void init(const std::vector<Point> &points, double radius, Vector3D bound_min, Vector3D bound_max);
    PivotTriangle find_seed_triangle();
    PivotTriangle pivot(PivotTriangle pt);
    std::vector<Point*> all_points;
    std::unordered_set<Point*> used;

    static double dist(const Point &p);
    int get_active_edge();
    PivotTriangle retrieve_active_edge(int index);
    void join(PivotTriangle e, Point* sigma_k, Point* new_center, int index);
    void glue(PivotTriangle ik);
    bool on_front(Point *k);
    bool not_used(Point *k);
    bool front_contains_edge(PivotTriangle t);
    void insert_edge(vector<PivotTriangle> edge);
    void mark_as_boundary(PivotTriangle e);

  private:
    // cell info for hashing
    double cell_width;

    struct CellIndex {
      unsigned long long int x_ind;
      unsigned long long int y_ind;
      unsigned long long int z_ind;

      // constructors
      CellIndex(unsigned long long int x, unsigned long long int y, unsigned long long int z) : x_ind( x ), y_ind( y ), z_ind( z ) { }

      // default constructor; be careful with this
      CellIndex() : x_ind( 0 ), y_ind( 0 ), z_ind( 0 ) { }

      // equality overrides
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

    std::vector<std::vector<PivotTriangle> > front;

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

    void create_spatial_grid(const std::vector<Point> &points);
    std::vector<Point *> neighborhood(double r, const Point &p);
    Vector3D circumcenter(const Point &a, const Point &b, const Point &c);
    bool valid_vertices(const Point &a, const Point &b, const Point &c);
    // overload ball_center to save some computation in the case that we
    // already know the triangle normal
    Point* ball_center(const Point &a, const Point &b, const Point &c);
    Point* ball_center(const Point &a, const Point &b, const Point &c, const Vector3D &normal);
    Vector3D naive_plane_normal(const Point &a, const Point &b, const Point &c);
    Vector3D correct_plane_normal(const Point &a, const Point &b, const Point &c);
    double ball_intersection(const Point &tc, double tr, const Point &ts, const Point &x);
    double angle_between(const Point &tc, const Point &ts, const Vector3D &i);
    int hash_position(const Point &p);
    int hash_cell(const CellIndex &c);
    CellIndex get_cell(const Point &p);
    Point* get_seed_candidate(const CellIndex &c);
    void increment_seed_cell();
    void calculate_normals();
    bool contains_edge(std::vector<PivotTriangle> vec, PivotTriangle e);
    void insert_edge(PivotTriangle e, VertexSpecifier v1, VertexSpecifier v2);
};

#endif //BALLPIVOT_H
