//
// Created by Rene Lee on 2019-04-23.
//

#include "CGL/CGL.h"
#include "ballPivot.h"
#include "point.h"

using namespace std;
using namespace CGL;

bool compare(Point *a, Point *b) {
  /* Sort in order of descending distance from *point
   * so that we can modify the list by popping from the back
   */
  BallPivot B;
  float dista = B.dist(*a);
  float distb = B.dist(*b);
  return dista > distb;
}

void BallPivot::init(std::vector <Point> points, float radius) {
  //this->used;
  this->unused = points;
  this->radius = radius;
  BallPivot::create_spatial_grid();
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
  while (!found_valid_triangle) {
    point = &unused[index];

    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    float h = hash_position(*point);
    if (map.find(h) != map.end()) {
      // TODO obtain a list of points in a (2 * rho)-neighborhood of *point,
      // or on the boundary of said neighborhood
      // (currently this just gets points in the same spatial partition)
      vector<Point *> *lst = map.at(h);

      // now, build potential seed triangles
      // organize lst in order of distance from point
      // such that closer points are at the back
      // TODO ask Rene about the lst->begin()+4
      // Rene: lst->begin()+4 was from an example I saw online
      std::sort(lst->begin(), lst->end(), compare);
//      std::sort(lst->begin()+4, lst->end(),
//              [x] (const Point& lhs, const Point& rhs) {
//                  return distance(x, lhs) < distance(x, rhs);
//              }
//              );


      // TODO: what happened to the hashmap of normals per vertex?
      Vector3D curr_normal = point->normal;
      // Stop when a valid seed triangle is found
      int i = 0;
      while (!found_valid_triangle && lst->size() >= 2) {

        // check that the triangle normal is consistent with the vertex normals, i.e. pointing outward
        // basically we need to check that all three vertices are pointing to the same side of the plane that the triangle creates
        Point *pointa = lst->at(i);
        Point *pointb = lst->at(i + 1);


        bool consistent = false;
        // 1. get the plane equation of the 3 points
        // 2. get the normal of the plane
        // 3. check that all 3 points are on the same side of the plane
        Vector3D CA = pointa->pos - point->pos;
        Vector3D CB = pointb->pos - point->pos;
        Vector3D plane_normal = cross(CA, CB);
        plane_normal.normalize();
        float k = -1 * dot(plane_normal, point->pos);

        // A, B, and C are points at the end of the normals corresponding to pointa, pointb, and point respectively
        Vector3D A = pointa->pos + pointa->normal;
        Vector3D B = pointb->pos + pointb->normal;
        Vector3D C = point->pos + point->normal;

        float sideA = dot(plane_normal, A) + k;
        float sideB = dot(plane_normal, B) + k;
        float sideC = dot(plane_normal, C) + k;

        if ((sideA > 0 && sideB > 0 && sideC > 0) || (sideA < 0 && sideB < 0 && sideC < 0)) {
          consistent = true;
        }

        // test that the p-ball with center in the outward half space touches all three vertices and contains no other data point
      }
    }

    index++;
    if (index >= unused.size() && !found_valid_triangle) {
      // No seed triangle was found!!
      std::vector<Point> empty;
      return empty;
    } else if (found_valid_triangle) {
      return triangle;
    }
  }
  return triangle;
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
   */

  Vector3D proj_center = circumcenter(a, b, c);
  double perp_dist = sqrt(pow(rho, 2) - (proj_center - c.pos).norm2());
  Vector3D plane = cross((a.pos - c.pos), (b.pos - c.pos));
  plane.normalize(); // Note: normalize() is a void function;
  Vector3D plane_normal = plane;

  // Ensure that the plane normal is pointing in the same direction as the
  // triangle normal (such that th calculated center of the sphere will be
  // in the correct half-space, "on top" of the triangle). All triangle normals
  // should be pointing to the same half-space, so we can simply check against
  // one of the vertices.
  if (dot(plane_normal, a.normal) < 0) {
    plane_normal *= -1;
  }

  return proj_center + perp_dist * plane_normal;
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

float BallPivot::distance(const Point &a, const Point &b) {
  Vector3D ab = a.pos - b.pos;
  return ab.norm();
}

float BallPivot::dist(const Point &a) {
  Vector3D ab = a.pos - point->pos;
  return ab.norm();
}


