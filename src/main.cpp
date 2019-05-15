#include "CGL/CGL.h"

#include "collada.h"
#include "meshEdit.h"
#include "bezierPatch.h"
#include "bezierCurve.h"
#include "mergeVertices.h"
#include "shaderUtils.h"
#include "rply.h"
#include "rplyfile.h"
#include "ballPivot.h"

#include <iostream>

using namespace std;
using namespace CGL;

#define msg(s) cerr << "[Collada Viewer] " << s << endl;

vector <Vector3D> vertices;
Vector3D vertex;

static int vertex_cb(p_ply_argument argument) {
  long eol;
  ply_get_argument_user_data(argument, NULL, &eol);
  double value = ply_get_argument_value(argument);
  if (eol == 0) {
    vertex = Vector3D(value, 0, 0);
  } else if (eol == 1) {
    vertex.y = value;
  } else {
    vertex.z = value;
    vertices.push_back(vertex);
  }
  return 1;
}

int loadFile(MeshEdit *collada_viewer, const char *path) {
  bool verbose = true;
  Scene *scene = new Scene();

  std::string path_str = path;
  if (path_str.substr(path_str.length() - 4, 4) == ".ply") {
    if (verbose) cout << "(main) Parsing ply file for points..." << flush;
    p_ply ply = ply_open(path, NULL, 0, NULL);
    p_ply_element element = NULL;
    int success = ply_read_header(ply);
    long nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 2);
    if (!ply_read(ply)) return 1;
    ply_close(ply);

    vector <Point> points;

    // track bounding box for spatial hashing
    double min_x = INF_D;
    double max_x = -INF_D;
    double min_y = INF_D;
    double max_y = -INF_D;
    double min_z = INF_D;
    double max_z = -INF_D;

    for (auto const &v : vertices) {
      if (v.x < min_x) {
        min_x = v.x;
      } else if (v.x > max_x) {
        max_x = v.x;
      }

      if (v.y < min_y) {
        min_y = v.y;
      } else if (v.y > max_y) {
        max_y = v.y;
      }

      if (v.z < min_z) {
        min_z = v.z;
      } else if (v.z > max_z) {
        max_z = v.z;
      }

      Point p = Point(v);
      points.push_back(p);
    }
    if (verbose) cout << " Done\n";

    double radius = 50;

    Vector3D bound_min = Vector3D(min_x, min_y, min_z);
    Vector3D bound_max = Vector3D(max_x, max_y, max_z);

    if (verbose) cout << "\n(main) min x: " << bound_min.x << flush;
    if (verbose) cout << "\n(main) max x: " << bound_max.x << flush;

    vector <BallPivot::PivotTriangle> triangles;

    // TODO write main loops for ball pivoting and output
    BallPivot pivot = BallPivot();
    pivot.init(points, radius, bound_min, bound_max);
    int index;
    while (true) {
      if (verbose) cout << "\n-----------------------------" << flush;
      index = pivot.get_active_edge();
      while (index != -1) {
        BallPivot::PivotTriangle *t = pivot.retrieve_active_edge(index);
        if (verbose) cout << "\n(main) Found active edge" << flush;

        BallPivot::PivotTriangle t_k = pivot.pivot(*t);
        Point *k = t_k.sigma_o;

        if (!t_k.empty && (pivot.not_used(*k) || pivot.on_front(*k))) {
          if (verbose) cout << "\n(main) Valid triangle found by pivoting" << flush;
          triangles.push_back(t_k);
          pivot.join(*t, k, t_k.center, index);
          cout << "\n(main) Joined" << flush;
          BallPivot::PivotTriangle ki = BallPivot::PivotTriangle(k, t->sigma_i, t->sigma_j, t_k.center);
          BallPivot::PivotTriangle jk = BallPivot::PivotTriangle(t->sigma_j, k, t->sigma_i, t_k.center);
          if (pivot.front_contains_edge(ki)) {
            pivot.glue(ki);
            if (verbose) cout << "\n(main) Glued" << flush;
          }
          if (pivot.front_contains_edge(jk)) {
            pivot.glue(jk);
            if (verbose) cout << "\n(main) Glued" << flush;
          }
        } else {
            pivot.mark_as_boundary(t);
            if (verbose) cout << "\n(main) Found Boundary Edge" << flush;
        }
        index = pivot.get_active_edge();
        if (verbose) cout << "\n(main) New active edge found" << flush;
      }

      if (verbose) cout << "\n(main) Calling seed triangle... " << flush;
      BallPivot::PivotTriangle seed_triangle = pivot.find_seed_triangle();
      if (verbose) cout << "\n(main) Done calling seed triangle\n" << flush;

      if (!seed_triangle.empty) {
        // output triangle
        triangles.push_back(seed_triangle);
        vector <BallPivot::PivotTriangle> edge_ij;
        vector <BallPivot::PivotTriangle> edge_jk;
        vector <BallPivot::PivotTriangle> edge_ki;
        BallPivot::PivotTriangle jk = BallPivot::PivotTriangle(seed_triangle.sigma_j, seed_triangle.sigma_o,
                                                               seed_triangle.sigma_i, seed_triangle.center);
        BallPivot::PivotTriangle ki = BallPivot::PivotTriangle(seed_triangle.sigma_o, seed_triangle.sigma_i,
                                                               seed_triangle.sigma_j, seed_triangle.center);
        edge_ij.push_back(seed_triangle);
        edge_jk.push_back(jk);
        edge_ki.push_back(ki);
        pivot.insert_edge(edge_ij);
        pivot.insert_edge(edge_jk);
        pivot.insert_edge(edge_ki);
        if (verbose) cout << "\n(main) Found seed triangle\n" << flush;
      } else {
        if (verbose) cout << "\n(main) Did not find seed triangle\n" << flush;
        break;
      }
    }

//    cout << "Hello" << flush;

    Camera *cam = new Camera();
    cam->type = CAMERA;
    Node ply_node;
    ply_node.instance = cam;
    scene->nodes.push_back(ply_node);

    Polymesh *mesh = new Polymesh();
    mesh->type = POLYMESH;
    ply_node.instance = mesh;
    scene->points = pivot.all_points;
    scene->nodes.push_back(ply_node);
    scene->triangles = triangles;
  } else if (path_str.substr(path_str.length() - 4, 4) == ".dae") {
    if (ColladaParser::load(path, scene) < 0) {
      delete scene;
      return -1;
    }
  } else if (path_str.substr(path_str.length() - 4, 4) == ".bez") {
    Camera *cam = new Camera();
    cam->type = CAMERA;
    Node node;
    node.instance = cam;
    scene->nodes.push_back(node);
    Polymesh *mesh = new Polymesh();

    FILE *file = fopen(path, "r");
    int n = 0;
    fscanf(file, "%d", &n);
    for (int i = 0; i < n; i++) {
      BezierPatch patch;
      patch.loadControlPoints(file);
      patch.add2mesh(mesh);
      mergeVertices(mesh);
    }
    fclose(file);

    mesh->type = POLYMESH;
    node.instance = mesh;
    scene->nodes.push_back(node);
  } else {
    return -1;
  }

  collada_viewer->load(scene);

  GLuint tex = makeTex("envmap/envmap.png");
  if (!tex) tex = makeTex("../envmap/envmap.png");
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, tex);
  glActiveTexture(GL_TEXTURE2);

  return 0;
}

int main(int argc, char **argv) {

  const char *path = argv[1];
  std::string path_str = path;

  //////////////////////////////
  // Bezier curve viewer code //
  //////////////////////////////

  if (path_str.substr(path_str.length() - 4, 4) == ".bzc") {
    // Each file contains a single Bezier curve's control points
    FILE *file = fopen(path, "r");

    int numControlPoints;
    fscanf(file, "%d", &numControlPoints);

    BezierCurve curve(numControlPoints);
    curve.loadControlPoints(file);
    fclose(file);

    // Create viewer
    Viewer viewer = Viewer();
    viewer.set_renderer(&curve);
    viewer.init();
    viewer.start();

    exit(EXIT_SUCCESS);

    return 0;
  }

  // create viewer
  Viewer viewer = Viewer();

  // create collada_viewer
  MeshEdit *collada_viewer = new MeshEdit();

  // set collada_viewer as renderer
  viewer.set_renderer(collada_viewer);

  // init viewer
  viewer.init();

  // load tests
  if (argc == 2) {
    if (loadFile(collada_viewer, argv[1]) < 0) exit(0);
  } else {
    msg("Usage: ./meshedit <path to scene file>");
    exit(0);
  }

  // start viewer
  viewer.start();

  return 0;
}
