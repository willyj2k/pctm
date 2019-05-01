//
// Created by Rene Lee on 2019-04-23.
//

#include "seedSelection.h"

void seedSelection::init(std::vector <Vector3D> points, std::map <Vector3D, Vector3D> normals, float radius) {
  this.used;
  this.unused = points;
  this.normals = normals;
  this.radius = radius;
  seedSelection::create_spatial_grid();
}

void seedSelection::create_spatial_grid() {
  for (const auto &entry : map) {
    delete (entry.second);
  }
  map.clear();

  for (int i = 0; i < unused.size(); i++) {
    Vector3D *p = &unused[i];
    float h = hash_position(p);
    if (map.find(h) == map.end()) {
      // does not exist already
      vector < Vector3D * > *lst = new vector<Vector3D *>();
      lst->emplace_back(p);
      map.insert({h, lst});
    } else {
      // exists
      vector < Vector3D * > *lst = map.at(h);
      lst->emplace_back(p);
      map.insert({h, lst});
    }
  }
}

std::vector <Vector3D> seedSelection::find_seed_triangle() {
  bool found_valid_triangle = false;
  // pick a point SIGMA that has not been used by the reconstructed triangulation;
  int index = 0;
  while (!found_valid_triangle) {
    point = &unused[index];

    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    float h = hash_position(*point);
    if (map.find(h) != map.end()) {
      vector < Vector3D * > *lst = map.at(h);
      // now we have a list of Vector3D that are in the same hashed 3D box
      // build potential seed triangles

      // organize lst in order of distance from point
      // use helper function
      std::sort (lst.begin()+4, lst.end(), compare);

      // check that the triangle normal is consistent with the vertex normals, i.e. pointing outward
      // basically we need to check that all three vertices are pointing to the same side of the plane that the triangle creates



      // test that the p-ball with center in the outward half space touches all three vertices and contains no other data point

      // Stop when a valid seed triangle is found
    }

    index++;
    if (index >= unused.size() && !found_valid_triangle) {
      // No seed triangle was found!!
      return nullptr;
    }
  }
}

Vector3D seedSelection::circumcenter(const Vector3D &a, const Vector3D &b, const Vector3D &c) {
  /* Returns the Cartesian coordinates of the circumcenter of the triangle
   * with vertices a, b, c
   */

  // obtain the barycentric coordinates of the circumcenter; formula from
  // https://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates
  double a2 = (c - b).norm2();
  double b2 = (c - a).norm2();
  double c2 = (b - a).norm2();
  double bary_a = a2 * (b2 + c2 - a2);
  double bary_b = b2 * (c2 + a2 - b2);
  double bary_c = c2 * (a2 + b2 - c2);

  // normalize barycentric coordinates so that they sum to 1
  double bary_sum = bary_a + bary_b + bary_c;
  bary_a /= bary_sum;
  bary_b /= bary_sum;
  bary_c /= bary_sum;

  // return Cartesian coordinates
  return bary_a * a + bary_b * b + bary_c * c;
}

Vector3D seedSelection::rho_center(double rho, const Vector3D &a, const Vector3D &b, const Vector3D &c) {
  /* Returns the Cartesian coordinates of the center of a sphere with radius 
   * rho that intersects the points a, b, c
   */

  Vector3D c = circumcenter(a, b, c);
  vector3D plane_normal = cross((a - c), (b - c)).normalize();
  // TODO ensure that plane normal points in the correct direction
  return c + rho * plane_normal;
}

float seedSelection::hash_position(Vector3D pos) {
  double w = 3 * width / (2 * radius);
  double h = 3 * height / (2 * radius);
  double t = max(w, h);
  // truncate the position to a specific 3D box
  double xpos = floor(pos[0] / w);
  double ypos = floor(pos[1] / h);
  double zpos = floor(pos[2] / t);

  return pow(113, 1) * xpos + pow(113, 2) * ypos + pow(113, 3) * zpos;
}

float seedSelection::distance(const Vector3D &a, const Vector3D &b) {
  Vector3D ab = a - b;
  return ab.norm();
}

bool seedSelection::compare(Vector3D *a, Vector3D *b) {
  float dista = distance(*a, *point);
  float distb = distance(*b, *point);
  return dista < distb;
}
