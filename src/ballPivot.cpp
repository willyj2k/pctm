//
// Created by Rene Lee on 2019-04-23.
//

#include "CGL/CGL.h"
#include "ballPivot.h"
#include "point.h"
#include <iostream>

using namespace std;
using namespace CGL;

Point *point;

bool compare(Point *a, Point *b) {
  /* Sort in order of descending distance from *point
   * so that we can modify the list by popping from the back
   */
  // TODO: we can get rid of "BallPivot B;" if we make the "dist" function static
  // TODO: this means "point" has to be static too
  // TODO: this might mess with our use of "point" elsewhere!
  float dista = BallPivot::dist(*a);
  float distb = BallPivot::dist(*b);
  return dista < distb;
}

void BallPivot::init(std::vector <Point> points, float radius) {
  //this->used;
  this->unused = points;
  this->radius = radius;
  BallPivot::create_spatial_grid();
  BallPivot::calculate_normals();
}

void BallPivot::create_spatial_grid() {
  for (const auto &entry : map) {
    delete (entry.second);
  }
  map.clear();

  for (int i = 0; i < unused.size(); i++) {
    Point *p = &unused[i];
    float h = hash_position(*p);
    if (map.find(h) == map.end()) {
      // does not exist already
      vector<Point *> *lst = new vector<Point *>();
      lst->push_back(p);
      //map.insert({h, lst});
      map.insert(std::make_pair(h, lst));
    } else {
      // exists
      vector<Point *> *lst = map.at(h);
      lst->push_back(p);
      //map.insert({h, lst});
      map.insert(std::make_pair(h, lst));
    }
  }
}

std::vector<Point> BallPivot::find_seed_triangle() {
  bool found_valid_triangle = false;
  bool consistent_normals;
  // pick a point SIGMA that has not been used by the reconstructed triangulation;
  int index = 0;
  std::vector<Point> triangle;
  while (!found_valid_triangle && index < unused.size()) {
    sigma = &unused[index];
    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    float h = hash_position(*sigma);

    if (map.find(h) != map.end()) {
      // TODO obtain a list of points in a (2 * rho)-neighborhood of *point,
      // or on the boundary of said neighborhood
      // (currently this just gets points in the same spatial partition)
      vector<Point *> *lst = map.at(h);

      // now, build potential seed triangles
      // organize lst in order of distance from point
      // such that closer points are at the back
      // std::sort(lst->begin(), lst->end(), &compare);

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
          // test that the p-ball with center in the outward half space touches
          // all three vertices and contains no other data point
          // TODO
        }
      }
    }
    ++index;
  }
  if (found_valid_triangle) {
    return triangle
  } else {
    // No seed triangle was found!!
    std::vector<Point> empty;
    return empty;
  }
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

Vector3D BallPivot::rho_center(double rho, const Point &a, const Point &b, const Point &c) {
  /* Returns the Cartesian coordinates of the center of a sphere with radius
   * rho that intersects the points a, b, c
   *
   * Assumes that a, b and c form a valid triangle.
   */
  // compute the projection of the sphere center onto the triangle abc
  Vector3D proj_center = circumcenter(a, b, c);

  // a, b and c lie on the surface of the sphere, so we can apply the
  // the Pythagorean theorem to find the perpendicular distance from the center
  // of the sphere to the triangle: (circumcenter - c)^2 + perp_dist^2 = rho^2
  double perp_dist = sqrt(pow(rho, 2) - (proj_center - c.pos).norm2());

  Vector3D plane_normal = correct_plane_normal(a, b, c);

  return proj_center + perp_dist * plane_normal;
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

float BallPivot::hash_position(const Point &p) {
  double w = 3 * width / (2 * radius);
  double h = 3 * height / (2 * radius);
  double t = max(w, h);
  // truncate the position to a specific 3D box
  double xpos = floor(p.pos.x / w);
  double ypos = floor(p.pos.y / h);
  double zpos = floor(p.pos.z / t);

  return pow(113, 1) * xpos + pow(113, 2) * ypos + pow(113, 3) * zpos;
}

void BallPivot::calculate_normals() {
    for (auto pair : map) {
        vector<Point *>* points = pair.second;
        Vector3D centroid;
        for (int curr = 0; curr < points->size(); curr++) {
            centroid = Vector3D();
            for (int i = 0; i < points->size(); i++) {
                if (i == curr) {
                    continue;
                }
                centroid += ((*points)[i])->pos;
            }
            if (points->size() > 1) {
                centroid = centroid / (points->size() - 1);
            }
            Vector3D mag = ((*points)[curr])->pos - centroid;
            Vector3D dir = mag.unit();
            ((*points)[curr])->normal = dir;
        }
    }
}

float BallPivot::distance(const Point &a, const Point &b) {
  Vector3D ab = a.pos - b.pos;
  return ab.norm();
}

float BallPivot::dist(const Point &a) {
  Vector3D ab = a.pos - point->pos;
  return ab.norm();
}


