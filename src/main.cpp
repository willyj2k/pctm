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

vector<Vector3D> vertices;
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

int loadFile(MeshEdit* collada_viewer, const char* path) {

  Scene* scene = new Scene();

  std::string path_str = path;
  if (path_str.substr(path_str.length()-4, 4) == ".ply")
  {
    cout << "Parsing ply file for points..." << flush;
    p_ply ply = ply_open(path, NULL, 0, NULL);
    p_ply_element element = NULL;
    int success = ply_read_header(ply);
    long nvertices = ply_set_read_cb(ply, "vertex", "x", vertex_cb, NULL, 0);
    ply_set_read_cb(ply, "vertex", "y", vertex_cb, NULL, 1);
    ply_set_read_cb(ply, "vertex", "z", vertex_cb, NULL, 2);
    if (!ply_read(ply)) return 1;
    ply_close(ply);

    vector<Point> points;
    for (int i = 0; i < vertices.size(); ++i) {
      Point p = Point(vertices[i], Vector3D(0, 0, 0));
      points.push_back(p);
    }
    cout << " Done\n";

    BallPivot pivot = BallPivot();
    pivot.init(points, 0.001);
    pivot.find_seed_triangle();

    Camera* cam = new Camera();
    cam->type = CAMERA;
    Node ply_node;
    ply_node.instance = cam;
    scene->nodes.push_back(ply_node);

    Polymesh* mesh = new Polymesh();
    mesh->type = POLYMESH;
    ply_node.instance = mesh;
    scene->points = vertices;
    scene->nodes.push_back(ply_node);
  }
  else if (path_str.substr(path_str.length()-4, 4) == ".dae")
  {
    if (ColladaParser::load(path, scene) < 0) {
      delete scene;
      return -1;
    }
  }
  else if (path_str.substr(path_str.length()-4, 4) == ".bez")
  {
    Camera* cam = new Camera();
    cam->type = CAMERA;
    Node node;
    node.instance = cam;
    scene->nodes.push_back(node);
    Polymesh* mesh = new Polymesh();

    FILE* file = fopen(path, "r");
    int n = 0;
    fscanf(file, "%d", &n);
    for (int i = 0; i < n; i++)
    {
      BezierPatch patch;
      patch.loadControlPoints(file);
      patch.add2mesh(mesh);
      mergeVertices(mesh);
    }
    fclose(file);

    mesh->type = POLYMESH;
    node.instance = mesh;
    scene->nodes.push_back(node);
  }
  else
  {
    return -1;
  }

  collada_viewer->load( scene );

  GLuint tex = makeTex("envmap/envmap.png");
  if(!tex) tex = makeTex("../envmap/envmap.png");
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, tex);
  glActiveTexture(GL_TEXTURE2);

  return 0;
}

int main( int argc, char** argv ) {

  const char* path = argv[1];
  std::string path_str = path;

  //////////////////////////////
  // Bezier curve viewer code //
  //////////////////////////////

  if (path_str.substr(path_str.length()-4, 4) == ".bzc")
  {
    // Each file contains a single Bezier curve's control points
    FILE* file = fopen(path, "r");

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
  MeshEdit* collada_viewer = new MeshEdit();

  // set collada_viewer as renderer
  viewer.set_renderer(collada_viewer);

  // init viewer
  viewer.init();

  // load tests
  if( argc == 2 ) {
    if (loadFile(collada_viewer, argv[1]) < 0) exit(0);
  } else {
    msg("Usage: ./meshedit <path to scene file>"); exit(0);
  }

  // start viewer
  viewer.start();

  return 0;
}
