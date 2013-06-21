#pragma once

#include <vector>
#include "Vector2D.h"

class PolygonalPath
{
private:
    std::vector<std::pair<Vector2D, float> > points; //space x time
    std::vector<Vector2D> tangents;
public:
    PolygonalPath();
    PolygonalPath(std::vector<std::pair<Vector2D,float> > points);
    PolygonalPath(std::vector<std::pair<Vector2D,float> > points,
                  const std::vector<Vector2D> &tangents);

    unsigned numberOfPoints() const { return points.size(); }

    inline std::pair<Vector2D, float> getPoint(unsigned index) const {
    	return points[index];
    }
    inline void addPoint(unsigned index, std::pair<Vector2D, float> newPoint, const Vector2D &tangent){
      points.insert((points.begin() + index), newPoint );
      tangents.insert((tangents.begin() + index), tangent);
    }
    inline void removerPoint(unsigned index) {
      points.erase(points.begin() + index);
    }
    const Vector2D &getTangent(unsigned index) const{
        return tangents.at(index);
    }

    std::string toString();
    float length();
};
