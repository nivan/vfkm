#include "PolygonalPath.h"

using namespace std;

PolygonalPath::PolygonalPath()
{
    this->points = std::vector<pair<Vector2D,float> >();
}

PolygonalPath::PolygonalPath(std::vector<pair<Vector2D,float> > points)
{
    this->points.assign(points.begin(), points.end());
    if (points.size() == 1)
        tangents.push_back(Vector2D(0,0));
    else {
        for (size_t i=1; i<points.size(); ++i) {
            tangents.push_back((points[i].first - points[i-1].first) * (1.0 / (points[i].second - points[i-1].second)));
        }
    }
}

PolygonalPath::PolygonalPath(std::vector<pair<Vector2D,float> > points,
                             const std::vector<Vector2D> &tangents_)
{
    this->points.assign(points.begin(),points.end());
    tangents = tangents_;
}

#include <sstream>

std::string PolygonalPath::toString(){
    stringstream ss;

    ss << "Curve ";

    int numberOfPoints = points.size();

    for(int i = 0 ; i < numberOfPoints ; ++i){
        pair<Vector2D,float> currentPoint = points.at(i);

        ss << "(" << currentPoint.first.toString() << ")," << currentPoint.second << ") ; ";
    }

    return ss.str();
}

float PolygonalPath::length(){
    float total = 0;
    int numberPoints = numberOfPoints();

    for(int i = 0 ; i < numberPoints - 1 ; ++i){
        Vector2D p0 = points.at(i).first;
        Vector2D p1 = points.at(i+1).first;

        p1.subtract(p0);

        total += p1.length();
    }

    return total;
}
