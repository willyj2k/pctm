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
  /* Sort in order of ascending distance from *sigma
   */
  // TODO: this means "sigma" has to be static too
  // TODO: this might mess with our use of "sigma" elsewhere!
  double dista = BallPivot::dist(*a);
  double distb = BallPivot::dist(*b);
  return dista < distb;
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
  this->unused.clear();
  cout << " Done\n";
  cout << "Calculating normals..." << flush;
  BallPivot::calculate_normals();
  this->all_points = this->unused;
  cout << " Done\n";
}

void BallPivot::create_spatial_grid() {
  for (const auto &entry : spatial_map) {
    delete (entry.second);
  }
  spatial_map.clear();

  for (int i = 0; i < unused.size(); i++) {
    Point *p = &unused[i];
    int h = hash_position(*p);
    if (spatial_map.find(h) == spatial_map.end()) {
      // does not already exist
      vector<Point *> *lst = new vector<Point *>();
      lst->push_back(p);
      spatial_map.insert(make_pair(h, lst));
    } else {
      // already exists
      vector<Point *> *lst = spatial_map.at(h);
      lst->push_back(p);
      spatial_map.insert(make_pair(h, lst));
    }
  }
}

vector<Point *> BallPivot::find_seed_triangle() {
  bool found_valid_triangle = false;
  bool consistent_normals;
  // pick a point SIGMA that has not been used by the reconstructed triangulation;
  CellIndex search_cell = CellIndex(0, 0, 0);
  int unused_index = 0;
  vector<Point *> triangle;
  while (!found_valid_triangle && unused_index < unused.size()) {
    // update
    sigma = &unused[unused_index];

    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    int h = hash_position(*sigma);

    if (spatial_map.find(h) != spatial_map.end()) {
      // obtain a list of points in a (2 * rho)-neighborhood of *point,
      // or on the boundary of said neighborhood
      // (currently this just gets points in the same spatial partition)
      vector<Point *> lst = neighborhood(2 * radius, *sigma);

      // now, build potential seed triangles
      // organize lst in order of distance from point
      // such that closer points are at the back
      sort(lst.begin(), lst.end(), compare);

      // Stop when a valid seed triangle is found
      for (int i = 1; !found_valid_triangle && i < lst.size(); ++i) {
        // check that the triangle normal is consistent with the vertex normals
        Point *sigma_a = lst.at(i - 1);
        Point *sigma_b = lst.at(i);

        // triangle_normal will be the zero vector if the points don't form a
        // valid triangle
        Vector3D triangle_normal = correct_plane_normal(*sigma, *sigma_a, *sigma_b);

        // test that there exists a p-ball with center in the outward half
        // space that touches all three vertices and contains no other data
        // point
        if (triangle_normal.norm2() > 0) {
          Point center = Point(ball_center(*sigma, *sigma_a, *sigma_b, triangle_normal));
          vector<Point *> r_neighborhood = neighborhood(radius, center);
          if (r_neighborhood.size() == 3) {
            // we don't neet to check membership in r_neighborhood because
            // sigma, sigma_a and sigma_b are already guaranteed to be distance
            // (radius) away from the center of the ball (assuming everything
            // is working correctly...)
            found_valid_triangle = true;
            triangle.push_back(sigma);
            triangle.push_back(sigma_a);
            triangle.push_back(sigma_b);
          }
        }
      }
    }
    ++unused_index;
  }
  if (found_valid_triangle) {
    return triangle;
  } else {
    // No seed triangle was found!!
    vector<Point *> empty;
    return empty;
  }
}

vector<Point *> BallPivot::neighborhood(double r, const Point &p) {
  /* Return a vector of pointers to points within an r-neighborhood of p */
  vector<Point *> r_neighborhood = vector<Point *>();
  int reach = ceil(r / cell_width);
  CellIndex c = get_cell(p);
  CellIndex cur_cell;
  int cur_hash;
  vector<Point *> *cur_points;

  // literally check all the cells that are possibly within reach...
  for (int x = (c.x_ind - reach); x <= (c.x_ind + reach); ++x) {
    for (int y = (c.y_ind - reach); y <= (c.y_ind + reach); ++y) {
      for (int z = (c.z_ind - reach); z <= (c.z_ind + reach); ++z) {
        cur_cell = CellIndex(x, y, z);
        cur_hash = hash_cell(cur_cell);
        if (spatial_map.find(cur_hash) != spatial_map.end()) {
          cur_points = spatial_map.at(cur_hash);
          for (auto const &q : *cur_points) {
            // we actually include the boundary of the neighborhood because
            // we're interested in consider spheres that could intersect
            // such points
            if ((q->pos - p.pos).norm() <= r) {
              r_neighborhood.push_back(q);
            }
          }
        }
      }
    }
  }
  return r_neighborhood;
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
  CellIndex cell = get_cell(p);
  return hash_cell(cell);
}

int BallPivot::hash_cell(const BallPivot::CellIndex &c) {
  return (c.x_ind + small_prime * (c.y_ind + small_prime * c.z_ind)) % large_prime;
}

BallPivot::CellIndex BallPivot::get_cell(const Point &p) {
  int x_ind = floor((p.pos.x - bound_min.x) / cell_width);
  int y_ind = floor((p.pos.y - bound_min.y) / cell_width);
  int z_ind = floor((p.pos.z - bound_min.z) / cell_width);
  return CellIndex(x_ind, y_ind, z_ind);
}

void BallPivot::calculate_normals() {

    // for (auto pair : map) {
    //     vector<Point *>* points = pair.second;
    //     Vector3D centroid;
    //     for (int curr = 0; curr < points->size(); curr++) {
    //         Vector3D pos = ((*points)[curr])->pos;
    //         centroid = Vector3D();
    //         for (int i = 0; i < points->size(); i++) {
    //             if (i == curr) {
    //                 continue;
    //             }
    //             centroid += ((*points)[i])->pos;
    //         }
    //         if (points->size() > 1) {
    //             centroid = centroid / (points->size() - 1);
    //         }
    //         Vector3D mag = pos - centroid;
    //         Vector3D dir = mag.unit();
    //         ((*points)[curr])->normal = mag;
    //         this->unused.push_back(((*points)[curr]));
    //     }
    // }
  // TODO verify that this rewrite works
  for (auto& pair : spatial_map) {
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

double BallPivot::dist(const Point &p) {
  Vector3D diff = p.pos - sigma->pos;
  return diff.norm();
}
