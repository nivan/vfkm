#ifndef GRID_H
#define GRID_H

#include <vector>
#include <set>
#include "Vector2D.h"
#include "PolygonalPath.h"
#include "Vector.h"

struct TriangularFace{
    int indexV1;
    int indexV2;
    int indexV3;
};

struct PointLocation {
    TriangularFace face;
    float barycentric_coords[3];
};

struct Segment {
    PointLocation endpoint[2];
    float time[2];
    int index;

    // remember that c . x actually returns Lambda . (v1, v2).
    // Lambda is the (constant) mixing matrix, so we keep it implicit in summand.

    inline void add_cx(Vector &summand, const Vector &field) const {
        float v1 =
            field[endpoint[0].face.indexV1] * endpoint[0].barycentric_coords[0] +
            field[endpoint[0].face.indexV2] * endpoint[0].barycentric_coords[1] +
            field[endpoint[0].face.indexV3] * endpoint[0].barycentric_coords[2];
        float v2 =
            field[endpoint[1].face.indexV1] * endpoint[1].barycentric_coords[0] +
            field[endpoint[1].face.indexV2] * endpoint[1].barycentric_coords[1] +
            field[endpoint[1].face.indexV3] * endpoint[1].barycentric_coords[2];
        summand[index]   += v1;
        summand[index+1] += v2;
    }
    /*inline void add_cx(float &v1, float &v2, const Vector &field) const {
        float v1 =
            field[endpoint[0].face.indexV1] * endpoint[0].barycentric_coords[0] +
            field[endpoint[0].face.indexV2] * endpoint[0].barycentric_coords[1] +
            field[endpoint[0].face.indexV3] * endpoint[0].barycentric_coords[2];
        float v2 =
            field[endpoint[1].face.indexV1] * endpoint[1].barycentric_coords[0] +
            field[endpoint[1].face.indexV2] * endpoint[1].barycentric_coords[1] +
            field[endpoint[1].face.indexV3] * endpoint[1].barycentric_coords[2];
        summand[index]   += v1;
        summand[index+1] += v2;
    }*/

    // remember that v is necessarily of the format Lambda . v
    inline void add_cTx(Vector &resulting_field, const Vector &v, float w=1.0f) const {
        // LT.L = [[1/3, 1/6], [1/6, 1/3]]
        // B LT.L = 
        float v1 = v[index], v2 = v[index+1];
        resulting_field[endpoint[0].face.indexV1] += w * v1 * endpoint[0].barycentric_coords[0] / 3.0;
        resulting_field[endpoint[0].face.indexV2] += w * v1 * endpoint[0].barycentric_coords[1] / 3.0;
        resulting_field[endpoint[0].face.indexV3] += w * v1 * endpoint[0].barycentric_coords[2] / 3.0;
        resulting_field[endpoint[0].face.indexV1] += w * v2 * endpoint[0].barycentric_coords[0] / 6.0;
        resulting_field[endpoint[0].face.indexV2] += w * v2 * endpoint[0].barycentric_coords[1] / 6.0;
        resulting_field[endpoint[0].face.indexV3] += w * v2 * endpoint[0].barycentric_coords[2] / 6.0;

        resulting_field[endpoint[1].face.indexV1] += w * v1 * endpoint[1].barycentric_coords[0] / 6.0;
        resulting_field[endpoint[1].face.indexV2] += w * v1 * endpoint[1].barycentric_coords[1] / 6.0;
        resulting_field[endpoint[1].face.indexV3] += w * v1 * endpoint[1].barycentric_coords[2] / 6.0;
        resulting_field[endpoint[1].face.indexV1] += w * v2 * endpoint[1].barycentric_coords[0] / 3.0;
        resulting_field[endpoint[1].face.indexV2] += w * v2 * endpoint[1].barycentric_coords[1] / 3.0;
        resulting_field[endpoint[1].face.indexV3] += w * v2 * endpoint[1].barycentric_coords[2] / 3.0;
    }


};

struct CurveDescription
{
    std::vector<Segment> segments;
    int index;
    float length;

    Vector rhsx, rhsy;

    inline void add_cTcx(Vector &resultX, const Vector &x, float k_global=1.0f) const {
        Vector v(2*segments.size());
        for(size_t j = 0; j < segments.size(); ++j) {
            float k = k_global * (segments[j].time[1] - segments[j].time[0]);
            segments[j].add_cx(v, x);
            segments[j].add_cTx(resultX, v, k);
        }
    }
};

typedef struct Intersection {
    PointLocation location;
    int indexV1;
    int indexV2;
    VECTOR_TYPE lambda; //interpolation factor in (1-alpha) v1 + alpha v2
    Vector2D dirAtIntersection;
    VECTOR_TYPE distanceBetweenCurveVertices;
} Intersection;

class Grid {
private:
    int m_resolutionX; // number of points on the sides of the grid
    int m_resolutionY; // number of points on the sides of the grid

    //domain
    float m_x;
    float m_y;
    float m_w;
    float m_h;

    float m_delta_x;
    float m_delta_y;

public:
    Grid(float x, float y, float w, float h, int resolutionX, int resolutionY);
    TriangularFace getFaceWherePointLies(const Vector2D &v) const;

    inline int getResolutionX() const { return m_resolutionX; }
    inline int getResolutionY() const { return m_resolutionY; }
    inline int getNumberOfFaces() const { return 2 * (m_resolutionX - 1) * (m_resolutionY - 1); }
    inline int getNumberOfFaces(int resolutionX, int resolutionY) const { return 2 * (resolutionX - 1) * (resolutionY - 1); }
    inline void getDomain(float &xD, float& yD,float &wD, float& hD) {
      xD = m_x; yD = m_y; wD = m_w; hD = m_h; }

    TriangularFace getFace(int index);
    TriangularFace getFace(int index, int resolutionX, int resolutionY);

    Vector2D toGrid(const Vector2D &world_point) const {
        return Vector2D((world_point.X() - m_x) / m_w * (m_resolutionX - 1.0),
                        (world_point.Y() - m_y) / m_h * (m_resolutionY - 1.0));                        
    }

    Vector2D toWorld(const Vector2D &grid_point) const {
        return Vector2D(grid_point.X() / (m_resolutionX - 1.0) * m_w + m_x,
                        grid_point.Y() / (m_resolutionY - 1.0) * m_h + m_y);
    }

    inline Vector2D getGridVertex(int index) const {
      return Vector2D(index % m_resolutionX, index / m_resolutionX);
    }
    inline Vector2D getVertex(int index) const {        
      return Vector2D(m_x + m_delta_x * (index % m_resolutionX), m_y + m_delta_y * (index / m_resolutionX));
    }
    inline Vector2D getVertex(int index,int resolutionX, int resolutionY) const {
        float delta_x = (float)m_w / (float)(resolutionX - 1);
        float delta_y = (float)m_h / (float)(resolutionY - 1);
      return Vector2D(m_x + delta_x * (index % resolutionX), m_y + delta_y * (index / resolutionX));
    }
    inline Vector2D getVertex(int rowIndex, int colIndex) const {
      return Vector2D(m_x + m_delta_x * colIndex, m_y + m_delta_y * rowIndex);
    }

    void computeVectorFieldImplicit(const Vector& vfXComponent,
                                    const Vector& vfYComponent,
                                    const Vector2D& point, Vector2D &result) const;

    // return barycentric_coords of point in face.
    PointLocation locate_point(const Vector2D &point) const;

    // If face is fixed, set only barycentric coordinates within face.
    void locate_point(PointLocation &face, const Vector2D &point) const;

    void multiplyByLaplacian(std::vector<Vector2D> &vectorField);

    void operator*(Vector&);
    void multiplyByLaplacian(Vector&, Vector&) const; //multiply both vectorfield by the laplacian
    void multiplyByLaplacian2(Vector &, Vector&); // Second vector stores the diagonal of L^T L for jacobi preconditioning.

    void clipLine(PolygonalPath &p) const;


    CurveDescription curve_description(const PolygonalPath &curve) const;

    void computeConstraints(const PolygonalPath& curve, std::vector<Intersection> & constraints) const;
    //should be called after calling clip planes
    //assumes that every segment of the line
    //lies in a face
    
    struct Inter {
        Vector2D grid_point;
        float u; // barycentric coordinate along the segment
        enum { Vertical, Horizontal, Diagonal, EndPoint } kind;
    };
    std::vector<Inter> clipAgainstVerticalLines  (const Inter &g1, const Inter &g2) const;
    std::vector<Inter> clipAgainstHorizontalLines(const Inter &g1, const Inter &g2) const;

    inline int vertexIndex(int x, int y) const {
        return y * m_resolutionX + x;
    }
    inline int faceIndex(int x, int y, bool is_bottom) const {
        int r = m_resolutionX - 1;
        return is_bottom ? (y * 2 * r + x) : (y * 2 * r + x + r);
    }
};

#endif // GRID_H
