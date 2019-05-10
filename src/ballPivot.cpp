//
// Created by Rene Lee on 2019-04-23.
//

#include "CGL/CGL.h"
#include "ballPivot.h"
#include "point.h"
#include <iostream>

using namespace CGL;
using std::vector;
using std::make_pair;
using std::sort;
using std::cout;
using std::flush;

static Point *sigma;

bool compare(Point *a, Point *b) {
  /* Sort in order of descending distance from *sigma
   * so that we can modify the list by popping from the back
   */
  // TODO: this means "sigma" has to be static too
  // TODO: this might mess with our use of "sigma" elsewhere!
  double dista = BallPivot::dist(*a);
  double distb = BallPivot::dist(*b);
  return dista > distb;
}

void BallPivot::init(vector<Point> points, double radius, Vector3D bound_min, Vector3D bound_max) {
  cout << "Initializing Ball Pivot member variables..." << flush;
  //this->used;
  this->unused = points;
  this->radius = radius;
  this->bound_min = bound_min;
  this->bound_max = bound_max;
  this->cell_width = 2 * radius;
  cout << " Done\n";
  
  cout << "Creating Spatial Grid..." << flush;
  BallPivot::create_spatial_grid();
  cout << " Done\n";

  cout << "Calculating vertex normals..." << flush;
  BallPivot::calculate_normals();
  cout << " Done\n";
}

void BallPivot::create_spatial_grid() {
  for (const auto &entry : map) {
    delete (entry.second);
  }
  map.clear();

  for (int i = 0; i < unused.size(); i++) {
    Point *p = &unused[i];
    double h = hash_position(*p);
    if (map.find(h) == map.end()) {
      // does not already exist
      vector<Point *> *lst = new vector<Point *>();
      lst->push_back(p);
      map.insert(make_pair(h, lst));
    } else {
      // already exists
      vector<Point *> *lst = map.at(h);
      lst->push_back(p);
      map.insert(make_pair(h, lst));
    }
  }
}

vector<Point> BallPivot::find_seed_triangle() {
  bool found_valid_triangle = false;
  bool consistent_normals;
  // pick a point SIGMA that has not been used by the reconstructed triangulation;
  int index = 0;
  vector<Point> triangle;
  while (!found_valid_triangle && index < unused.size()) {
    // update 
    sigma = &unused[index];

    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    double h = hash_position(*sigma);

    if (map.find(h) != map.end()) {
      // TODO obtain a list of points in a (2 * rho)-neighborhood of *point,
      // or on the boundary of said neighborhood
      // (currently this just gets points in the same spatial partition)
      vector<Point *> *lst = map.at(h);

      // now, build potential seed triangles
      // organize lst in order of distance from point
      // such that closer points are at the back
      sort(lst->begin(), lst->end(), compare);

      // Stop when a valid seed triangle is found
      int i = 0;
      while (!found_valid_triangle && lst->size() >= 2) {
        // check that the triangle normal is consistent with the vertex normals
        Point *sigma_a = lst->at(i);
        Point *sigma_b = lst->at(i + 1);

        // triangle_normal will be the zero vector if the points don't form a
        // valid triangle
        Vector3D triangle_normal = correct_plane_normal(*sigma, *sigma_a, *sigma_b);

        if (triangle_normal.norm2() > 0) {
          // TODO
          // test that there exists a p-ball with center in the outward half
          // space that touches all three vertices and contains no other data
          // point
          Vector3D center = ball_center(*sigma, *sigma_a, *sigma_b, triangle_normal);
        }
      }
    }
    ++index;
  }
  if (found_valid_triangle) {
    return triangle;
  } else {
    // No seed triangle was found!!
    vector<Point> empty;
    return empty;
  }
}

vector<Point *> BallPivot::neighborhood(double r, const Point &p) {
  /* Return a vector of pointers to points within an r-neighborhood of p */
  // TODO
  int reach = ceil(r / cell_width);
  return vector<Point *>();
}

Vector3D BallPivot::circumcenter(const Point &a, const Point &b, const Point &c) {
  /* Returns the Cartesian coordinates of the circumcenter of the triangle
   * with vertices a, b, c
   */

  // obtain the barycentric coordinates of the circumcenter; formula from
  // https://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates
  double a2 = (c.pos - b.pos).norm2();
  double b2 = (c.pos - a.pos).norm2();
  double c2 = (b.pos - a.pos).norm2();
  double bary_a = a2 * (b2 + c2 - a2);
  double bary_b = b2 * (c2 + a2 - b2);
  double bary_c = c2 * (a2 + b2 - c2);

  // normalize barycentric coordinates so that they sum to 1
  double bary_sum = bary_a + bary_b + bary_c;
  bary_a /= bary_sum;
  bary_b /= bary_sum;
  bary_c /= bary_sum;

  // return Cartesian coordinates
  return bary_a * a.pos + bary_b * b.pos + bary_c * c.pos;
}

Vector3D BallPivot::ball_center(const Point &a, const Point &b, const Point &c, const Vector3D &normal) {
  /* Returns the Cartesian coordinates of the center of a sphere with radius
   * this->radius that intersects the points a, b, c.
   *
   * Assumes that a, b and c form a valid triangle, and takes the triangle's
   * normal as an input (just to save some redundant calculation).
   */

  // compute the projection of the sphere center onto the triangle abc
  Vector3D proj_center = circumcenter(a, b, c);

  // a, b and c lie on the surface of the sphere, so we can apply the
  // the Pythagorean theorem to find the perpendicular distance from the center
  // of the sphere to the triangle: (circumcenter - c)^2 + perp_dist^2 = radius^2
  double perp_dist = sqrt(pow(radius, 2) - (proj_center - c.pos).norm2());

  return proj_center + perp_dist * normal;
}

Vector3D BallPivot::naive_plane_normal(const Point &a, const Point &b, const Point &c) {
  /* Returns the normal (unit) vector of the plane defined by a, b, and c
   *
   * NOTE 1: (IMPORTANT) This could fail if a, b and c are colinear, in which
   * case the cross product will be 0 (and they don't define a triangle anyway)
   *
   * On the other hand, we can use this to our advantage in order to check
   * whether the points define a valid triangle
   *
   * NOTE 2: The vector will not be oriented "correctly" in general; i.e., its
   * sign could be incorrect
   */
  Vector3D normal = cross((b.pos - a.pos), (c.pos - a.pos));
  return normal.unit();
}

Vector3D BallPivot::correct_plane_normal(const Point &a, const Point &b, const Point &c) {
  /* Returns the correctly oriented plane normal, or the zero vector if no such
   * vector exists
   *
   * Note: We can therefore use this to also verify that the vertex normals
   * (of a, b and c) are aligned.
   */
  Vector3D naive_normal = naive_plane_normal(a, b, c);

  if (naive_normal.norm2() == 0) {
    return naive_normal;
  }

  double a_dot = dot(a.normal, naive_normal);
  double b_dot = dot(b.normal, naive_normal);
  double c_dot = dot(c.normal, naive_normal);
  bool all_pos = a_dot >= 0 && b_dot >= 0 && c_dot >= 0;
  bool all_neg = a_dot >= 0 && b_dot <= 0 && c_dot <= 0;

  if (all_pos) {
    return naive_normal;
  } else if (all_neg) {
    return -1 * naive_normal;
  } else {
    return Vector3D(0, 0, 0);
  }
}

int BallPivot::hash_position(const Point &p) {
  // divide the bounding box in to cubic cells with side length 2 * radius
  // truncate the position of p to a specific 3D box
  cell_index cell = get_cell(p); 
  return (cell.x_ind + small_prime * (cell.y_ind + small_prime * cell.z_ind)) % large_prime;
}

BallPivot::cell_index BallPivot::get_cell(const Point &p) {
  int x_ind = floor((p.pos.x - bound_min.x) / cell_width);
  int y_ind = floor((p.pos.y - bound_min.y) / cell_width);
  int z_ind = floor((p.pos.z - bound_min.z) / cell_width);
  return cell_index(x_ind, y_ind, z_ind);
}

void BallPivot::calculate_normals() {
  // TODO verify that this rewrite works
  for (auto& pair : map) {
    vector<Point *> *points = pair.second;
    // first calculate the centroid of the cell (cube)
    // by taking the average of the position vectors
    Vector3D centroid = Vector3D(0, 0, 0);
    for (auto const &point : *points) {
      centroid += point->pos;
    }

    // now assign the normals of each point to be the difference
    // between the point's position and the average position of
    // the rest of the points in the cell
    Vector3D avg_other_pos;
    double num_other = (points->size() > 1) ? points->size() - 1 : 1;
    for (auto &point : *points) {
      avg_other_pos = (centroid - point->pos) / num_other;
      point->normal = (point->pos - avg_other_pos).unit();
    }
  }
}

double BallPivot::dist(const Point &p) {
  Vector3D diff = p.pos - sigma->pos;
  return diff.norm();
}

