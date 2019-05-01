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
    Vector3D
    p * = &unused[index];

    // consider all pairs of points in its neighborhood
    // first get the neighborhood, aka use spatial map
    float h = hash_position(*p);
    if (map.find(h) != map.end()) {
      vector < Vector3D * > *lst = map.at(h);
      // now we have a list of Vector3D that are in the same hashed 3D box
      // build potential seed triangles

      // organize lst in order of distance from p
      // use helper function


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


// This helper function is used to help make sure we aren't using the same point for a triangle
bool seedSelection::equal_positions(Vector3D pos1, Vector3D pos2) {
  return pos1[0] == pos2[0] && pos1[1] == pos2[1] && pos1[2] == pos2[2];
}

float seedSelection::distance(Vector3D a, Vector3D b) {
  Vector3D ab = a - b;
  return sqrt(pow(ab[0], 2) + pow(ab[1], 2) + pow(ab[2], 2));
}

bool seedSelection::compare(Vector3D a, Vector3D b) {
  
}