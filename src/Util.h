#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include "PolygonalPath.h"
#include "Vector.h"

typedef struct Cluster{
    std::string name;
    std::vector<Cluster*> children;
    std::vector<int> indices;
    std::vector<float> curveErrors;
    std::pair<Vector*,Vector*> vectorField;
    Cluster* parent;
    float error;
    float maxError;

    void clearChildren();
} Cluster;

class Util
{
public:
    Util();

    static void loadCurves(std::string filename, std::vector<PolygonalPath>&);
    static void loadCurves(std::string filename, std::vector<PolygonalPath>&,
                           float &xmin, float &xmax, float &ymin, float &ymax,
                           float &tmin, float &tmax);
    static void loadCurvesAndProject(std::string filename, std::vector<PolygonalPath>&,
                                     float &xmin, float &xmax, float &ymin, float &ymax,
                                     float &tmin, float &tmax);

    static void to_mercator(const float &lat, const float &lon, float &xMerc, float &yMerc);
};

#endif // UTIL_H
