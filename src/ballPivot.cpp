//
// Created by Rene Lee on 2019-04-23.
//

#include "CGL/CGL.h"
#include "ballPivot.h"
#include "point.h"
#include <iostream>
#include <unordered_set>
#include <cmath>
#include <algorithm>

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

void BallPivot::init(const vector<Point> &points, double radius, Vector3D bound_min, Vector3D bound_max) {
  cout << "Initializing Ball Pivot member variables..." << flush;
  this->radius = radius;
  this->bound_min = bound_min;
  this->bound_max = bound_max;
  this->seed_cell = CellIndex(0, 0, 0);
  this->max_cell = get_cell(Point(bound_max));
  // add 1 to max_cell_z because we'll traverse the cells such that z values
  // or contiguous
  this->max_cell.z_ind += 1;
  this->cell_width = 2 * radius;
  cout << " Done\n";

  cout << "Creating Spatial Grid..." << flush;
  BallPivot::create_spatial_grid(points);
  cout << " Done\n";
  cout << "Calculating normals..." << flush;
  BallPivot::calculate_normals();
  cout << " Done\n";
}

void BallPivot::create_spatial_grid(const vector<Point> &points) {
  for (const auto &entry : spatial_map) {
    delete (entry.second);
  }
  spatial_map.clear();

  for (Point p : points) {
    int h = hash_position(p);
    if (spatial_map.find(h) == spatial_map.end()) {
      // does not already exist
      vector<Point> *lst = new vector<Point>();
      lst->push_back(p);
      spatial_map.insert(make_pair(h, lst));
    } else {
      // already exists
      vector<Point> *lst = spatial_map.at(h);
      lst->push_back(p);
      spatial_map.insert(make_pair(h, lst));
    }
  }
}

BallPivot::PivotTriangle BallPivot::find_seed_triangle() {
  bool found_valid_triangle = false;
  PivotTriangle triangle;
  // pick a point SIGMA that has not been used by the reconstructed triangulation;
  while (!found_valid_triangle && seed_cell != max_cell) {
    int h = hash_cell(seed_cell);

    if (processed_cells.find(h) == processed_cells.end()) {
      sigma = get_seed_candidate(seed_cell);

      // consider all pairs of points in its neighborhood
      // first get the neighborhood, aka use spatial map
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

          if (valid_vertices(*sigma, *sigma_a, *sigma_b)) {
            // triangle_normal will be the zero vector if the points don't form a
            // valid triangle
            Vector3D triangle_normal = correct_plane_normal(*sigma, *sigma_a, *sigma_b);

            // test that there exists a p-ball with center in the outward half
            // space that touches all three vertices and contains no other data
            // point
            if (triangle_normal.norm2() > 0) {
              Point center = *ball_center(*sigma, *sigma_a, *sigma_b, triangle_normal);
              vector<Point *> r_neighborhood = neighborhood(radius, center);
              if (r_neighborhood.size() == 3) {
                // we don't neet to check membership in r_neighborhood because
                // sigma, sigma_a and sigma_b are already guaranteed to be distance
                // (radius) away from the center of the ball (assuming everything
                // is working correctly...)
                found_valid_triangle = true;
                triangle = PivotTriangle(sigma, sigma_a, sigma_b, &center);
                used.insert(sigma);
                used.insert(sigma_a);
                used.insert(sigma_b);
              }
            }
          }
        }
      }
      // put this here because apparently we only want to consider one
      // candidate *vertex* per cell, rather than one seed triangle per cell
      processed_cells.insert(h);
    }
    increment_seed_cell();
  }

  return triangle;
}

BallPivot::PivotTriangle BallPivot::pivot(BallPivot::PivotTriangle pt) {
  /* Return the triangle vertices of the new triangle and the center of the
   * corresponding ball found by pivoting if such a triangle is found.
   *
   * Always pivots around the first two edges in pt; pt.sigma_i and pt.sigma_j.
   *
   * Takes in the vertices and center of the corresponding ball for the
   * previous triangle.
   */
  if (pt.empty) {
    return pt;
  }
  Vector3D mid_ij = (pt.sigma_i->pos + pt.sigma_j->pos) / 2.0;
  Vector3D tri_normal = correct_plane_normal(*(pt.sigma_i), *(pt.sigma_j), *(pt.sigma_o));
  Vector3D proj_center = circumcenter(*(pt.sigma_i), *(pt.sigma_j), *(pt.sigma_o));
  Vector3D rotation_axis = cross(tri_normal, mid_ij - proj_center).unit();

  Point m = Point(mid_ij, rotation_axis);
  double trajectory_radius = (pt.center->pos - m.pos).norm();

  vector<Point *> candidates = neighborhood(2 * radius, m);
  Point *first_hit;
  Point *first_center;
  double min_theta = INF_D;
  for (Point *sigma_x : candidates) {
    if (valid_vertices(*(pt.sigma_i), *(pt.sigma_j), *sigma_x)) {
      Point *c_x = ball_center(*(pt.sigma_i), *(pt.sigma_j), *c_x);
      double theta = ball_intersection(m, trajectory_radius, *(pt.center), *c_x);
       if (theta > 0 && theta < 2 * PI && theta < min_theta) {
         min_theta = theta;
         first_hit = sigma_x;
         first_center = c_x;
       }
    }
  }
  if (first_hit != NULL) {
    return PivotTriangle(pt.sigma_i, pt.sigma_j, first_hit, first_center);
  } else {
    return PivotTriangle();
  }
}

double BallPivot::ball_intersection(const Point &tc, double tr, const Point &ts, const Point &x) {
  /* Returns the angle along the circular trajectory defined by center tc,
   * radius tr and starting point ts that the p-ball hits the point x. If the
   * ball never hits x, the this returns 0 instead.
   *
   * Assumes that the normal of tc is orthogonal to the trajectory plane and is
   * oriented so that the trajectory goes counterclockwise when facing against
   * the normal (i.e., if the normal is pointing "up" out of the face of
   * the clock).
   *
   * Implementation adapted from https://gamedev.stackexchange.com/questions/
   * 75756/sphere-sphere-intersection-and-circle-sphere-intersection
   */
  double d = abs(dot(tc.normal, (x.pos - tc.pos)));
  if (d > radius) {
    // trajectory plane doesn't intersect the ball
    return 0;
  }

  Vector3D intersection;

  // center of the ball centered at x, projected onto the trajectory plane
  Vector3D x_pc = x.pos + d * tc.normal;

  if (d == radius) {
    // the trajectory plane is tangent to the ball, so x_pc is the only
    // intersection point
    if ((x_pc - tc.pos).norm() == tr) {
      return angle_between(tc, ts, x_pc);
    } else {
      return 0;
    }
  }

  // radius of the circular slice of the ball in the trajectory plane
  double x_pr = sqrt(radius * radius  - d * d);
  
  // now we do circle-circle intersection
  double d_p = (tc.pos - x_pc).norm();
  if (tr + x_pr > d_p) {
    // circles too far away
    return 0;

  } else if (d_p + min(tr, x_pr) < max(tr, x_pr)) {
    // one circle entirely contained in the other
    return 0;

  } else if (tr + x_pr == d_p) {
    // circles are tangent (exterior)
    intersection = tc.pos + tr * (tc.pos - x_pc).unit();
    return angle_between(tc, ts, intersection);

  } else if (d_p + tr == x_pr) {
    // circles are tangent, trajectory inside ball
    intersection = x_pc + x_pr * (tc.pos - x_pc).unit();
    return angle_between(tc, ts, intersection);

  } else if (d_p + x_pr == tr) {
    // circles are tangent, ball inside trajectory
    intersection = tc.pos + tr * (x_pc - tc.pos).unit();
    return angle_between(tc, ts, intersection);

  } else {
    // circles intersect at two points
    double h = 0.5 + (tr * tr - x_pr * x_pr) / (2 * d_p * d_p);
    Vector3D c_i = tc.pos + h * (x_pc - tc.pos);
    double r_i = sqrt(tr * tr - h * h * d_p * d_p);
    Vector3D r_i_dir = cross(tc.normal, c_i).unit();
    Vector3D intersection1 = c_i + r_i * r_i_dir;
    Vector3D intersection2 = c_i - r_i * r_i_dir;
    double theta1 = angle_between(tc, ts, intersection1);
    double theta2 = angle_between(tc, ts, intersection2);
    return min(theta1, theta2);
  }
}

double BallPivot::angle_between(const Point &tc, const Point &ts, const Vector3D &i) {
  /* Returns the angle between (tc.pos - tc.pos) and (i - tc.pos); the angle is
   * measured counterclockwise from ts about tc (where the normal of tc is
   * facing out of the clock face)
   */
  Vector3D start = (ts.pos - tc.pos);
  Vector3D end = (i - tc.pos);
  return atan2(dot(cross(start, end), tc.normal), dot(start, end));
}

vector<Point *> BallPivot::neighborhood(double r, const Point &p) {
  /* Return a vector of pointers to points within an r-neighborhood of p */
  vector<Point *> r_neighborhood = vector<Point *>();
  int reach = ceil(r / cell_width);
  CellIndex c = get_cell(p);
  CellIndex cur_cell;
  int cur_hash;
  vector<Point> *cur_points;

  // literally check all the cells that are possibly within reach...
  for (int x = (c.x_ind - reach); x <= (c.x_ind + reach); ++x) {
    for (int y = (c.y_ind - reach); y <= (c.y_ind + reach); ++y) {
      for (int z = (c.z_ind - reach); z <= (c.z_ind + reach); ++z) {
        cur_cell = CellIndex(x, y, z);
        cur_hash = hash_cell(cur_cell);
        if (spatial_map.find(cur_hash) != spatial_map.end()) {
          cur_points = spatial_map.at(cur_hash);
          for (auto &q : *cur_points) {
            // we actually include the boundary of the neighborhood because
            // we're interested in consider spheres that could intersect
            // such points
            if ((q.pos - p.pos).norm() <= r) {
              r_neighborhood.push_back(&q);
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

bool BallPivot::valid_vertices(const Point &a, const Point &b, const Point &c) {
  /* Returns true if the input vertices can be touched by a ball of radius
   * this->radius.
   */
  Vector3D proj_center = circumcenter(a, b, c);
  double a_dist = (a.pos - proj_center).norm();
  if (a_dist > radius) return false;

  double b_dist = (b.pos - proj_center).norm();
  if (b_dist > radius) return false;

  double c_dist = (c.pos - proj_center).norm();
  if (c_dist > radius) return false;

  if (correct_plane_normal(a, b, c).norm2() == 0) return false;

  return true;
}

Point* BallPivot::ball_center(const Point &a, const Point &b, const Point &c) {
  const Vector3D normal = correct_plane_normal(a, b, c);
  return ball_center(a, b, c, normal);
}

Point* BallPivot::ball_center(const Point &a, const Point &b, const Point &c, const Vector3D &normal) {
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

  return new Point(proj_center + perp_dist * normal, normal);
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
   * vector exists.
   *
   * Note: We can therefore use this to verify that the vertex normals
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

Point* BallPivot::get_seed_candidate(const CellIndex &c) {
  /* Return the point in the cell aligned most strongly
   * with the average normal
   */
  int c_hash = hash_cell(c);
  if (spatial_map.find(c_hash) != spatial_map.end()) {
    // find average normal
    vector<Point> *points = spatial_map.at(c_hash);
    Vector3D avg_normal = Vector3D(0, 0, 0);
    for (const auto &p : *points) {
      avg_normal += p.normal;
    }
    avg_normal /= points->size();

    double max_dot = -INF_D;
    Point *best_candidate;
    for (Point &p : *points) {
      double cur_dot = dot(p.normal, avg_normal);
      if (cur_dot > max_dot) {
        max_dot = cur_dot;
        best_candidate = &p;
      }
    }
    return best_candidate;
  } else {
    return NULL;
  }
}

void BallPivot::increment_seed_cell() {
  if (seed_cell.z_ind >= max_cell.z_ind) {
    seed_cell.z_ind = 0;
    seed_cell.y_ind += 1;
    if (seed_cell.y_ind >= max_cell.y_ind) {
      seed_cell.y_ind = 0;
      seed_cell.x_ind += 1;
      if (seed_cell.x_ind >= max_cell.x_ind) {
        // at this point all cells should
        // have been traversed
      }
    } else {
      seed_cell.y_ind += 1;
    }
  } else {
    seed_cell.z_ind += 1;
  }
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
      for (const auto &pair : spatial_map) {
        vector<Point> *points = pair.second;
        // first calculate the centroid of the cell (cube)
        // by taking the average of the position vectors
        Vector3D centroid = Vector3D(0, 0, 0);
        for (const auto &point : *points) {
          centroid += point.pos;
        }

        // now assign the normals of each point to be the difference
        // between the point's position and the average position of
        // the rest of the points in the cell
        Vector3D avg_other_pos;
        double num_other = (points->size() > 1) ? points->size() - 1 : 1;
        for (auto &point : *points) {
          avg_other_pos = (centroid - point.pos) / num_other;
          point.normal = (point.pos - avg_other_pos);
          this->all_points.push_back(&point);
        }
    }
}

double BallPivot::dist(const Point &p) {
  Vector3D diff = p.pos - sigma->pos;
  return diff.norm();
}

bool compare_3D(Vector3D a, Vector3D b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

void BallPivot::join(PivotTriangle e, Point* sigma_k, Point* new_center, int index) {
    this->front[index].pop_back();
    PivotTriangle ik = PivotTriangle(e.sigma_i, sigma_k, e.sigma_j, new_center);
    PivotTriangle kj = PivotTriangle(sigma_k, e.sigma_j, e.sigma_i, new_center);
    this->front[index].push_back(ik);
    this->front[index].push_back(kj);
}

bool compare_edge(BallPivot::PivotTriangle e1, BallPivot::PivotTriangle e2) {
    return compare_3D(e1.sigma_i->pos, e2.sigma_i->pos) && compare_3D(e1.sigma_j->pos, e2.sigma_j->pos);
}

bool BallPivot::contains_edge(vector<PivotTriangle> vec, PivotTriangle e) {
    for (int i = 0; i < vec.size(); i++) {
        if (compare_edge(e, vec[i])) {
            return true;
        }
    }
    return false;
}

void BallPivot::glue(PivotTriangle ij) {
    PivotTriangle ji = PivotTriangle(ij.sigma_j, ij.sigma_i, ij.sigma_o, ij.center);
    int loop_index1 = 0;
    int loop_index2 = 0;
    int index1 = 0;
    int index2 = 0;

    for (int i = 0; i < front.size(); i++) {
        if (contains_edge(front[i], ij)) {
            loop_index1 = i;
            for (int j = 0; j < front[i].size(); j++) {
                if (compare_edge(front[i][j], ij)) {
                    index1 = j;
                }
            }
        }
        if (contains_edge(front[i], ji)) {
            loop_index2 = i;
            for (int j = 0; j < front[i].size(); j++) {
                if (compare_edge(front[i][j], ji)) {
                    index2 = j;
                }
            }
        }
    }

    //no duplicate edges contained- no need for gluing
    if (!(contains_edge(front[loop_index1], ij) && contains_edge(front[loop_index2], ji))) {
        return;
    }

    //edges belong to the same loop
    if (loop_index1 == loop_index2) {
        //edges form entirety of loop (Scenario a)
        if (front[loop_index1].size() == 2) {
            front.erase(front.begin() + (loop_index1));
        } else {
            //edges form a loop and are adjacent
            if (std::abs(index1 - index2) == 1) {
                front[loop_index1].erase(front[loop_index1].begin() + (index1));
                if (index1 < index2) {
                    front[loop_index1].erase(front[loop_index1].begin() + (index2 - 1));
                } else {
                    front[loop_index1].erase(front[loop_index1].begin() + (index2));
                }
            } else {
            //edges form a loop and are not adjacent
                vector<PivotTriangle> loop1;
                vector<PivotTriangle> loop2;
                if (index1 < index2) {
                    for (int i = index1 + 1; i < index2; i++) {
                        loop1.push_back(front[loop_index1][i]);
                    }
                    for (int i = index2 + 1; i < front[loop_index1].size(); i++) {
                        loop2.push_back(front[loop_index1][i]);
                    }
                    for (int i = 0; i < index1; i++) {
                        loop2.push_back(front[loop_index1][i]);
                    }
                } else {
                    for (int i = index2 + 1; i < index1; i++) {
                        loop1.push_back(front[loop_index1][i]);
                    }
                    for (int i = index1 + 1; i < front[loop_index1].size(); i++) {
                        loop2.push_back(front[loop_index1][i]);
                    }
                    for (int i = 0; i < index2; i++) {
                        loop2.push_back(front[loop_index1][i]);
                    }
                }
                front.erase(front.begin() + (loop_index1));
                front.push_back(loop1);
                front.push_back(loop2);
            }
        }
    } else {
        //edges are in different loops
        vector<PivotTriangle> loop1;
        for (int i = 0; i < index1; i++) {
            loop1.push_back(front[loop_index1][i]);
        }
        for (int i = index2 + 1; i < front[loop_index2].size(); i++) {
            loop1.push_back(front[loop_index2][i]);
        }
        for (int i = 0; i < index2; i++) {
            loop1.push_back(front[loop_index2][i]);
        }
        for (int i = index1 + 1; i < front[loop_index1].size(); i++) {
            loop1.push_back(front[loop_index1][i]);
        }
        front.erase(front.begin() + (loop_index1));
        if (loop_index1 < loop_index2) {
            front.erase(front.begin() + (loop_index2 - 1));
        } else {
            front.erase(front.begin() + (loop_index2));
        }
        front.push_back(loop1);
    }
}

bool BallPivot::on_front(Point k) {
    bool internal_mesh_vertex = false;
    for (int i = 0; i < front.size(); i++) {
        for (int j = 0; j < front.size(); j++) {
            if (compare_3D(front[i][j].sigma_i->pos, k.pos) || compare_3D(front[i][j].sigma_j->pos, k.pos)) {
                internal_mesh_vertex = true;
            }
        }
    }
    return internal_mesh_vertex;
}

bool BallPivot::not_used(Point k) {
    return !(used.find(&k) == used.end());
}

void mark_as_boundary(BallPivot::PivotTriangle e) {
    e.isBoundary = true;
}
