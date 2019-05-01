//
// Created by Rene Lee on 2019-04-23.
//

#ifndef COLLADAVIEWER_SEEDSELECTION_H
#define COLLADAVIEWER_SEEDSELECTION_H


class seedSelection {
  public:
    void init(std::vector<Point> points, float radius);
    std::vector<Point> find_seed_triangle();

  private:
    // std::vector of used points
    std::vector<Point> used;

    // std::vector of unused points
    std::vector<Point> unused;

    // spatial map
    unordered_map<float, vector<Point *> *> map;

    // float radius
    float radius;

    double width;
    double height;

    Point *point;

    void create_spatial_grid ();
    float hash_position(Vector3D pos);
    bool equal_positions(Vector3D pos1, Vector3D pos2);
    float distance(Vector3D a, Vector3D b);
};


#endif //COLLADAVIEWER_SEEDSELECTION_H
