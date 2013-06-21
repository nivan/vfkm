#ifndef OPTIMIZER_H
#define OPTIMIZER_H

//#include <vector>
#include "Grid.h"

struct ProblemSettings
{
    Grid &grid;
    const std::vector<int> &curveIndices;
    const std::vector<CurveDescription> &curve_descriptions;
    float totalCurveLength;
    float smoothnessWeight;
    ProblemSettings(Grid &g,
                    const std::vector<int> &i,
                    const std::vector<CurveDescription> &cd,
                    float tcl, float sw):
        grid(g),
        curveIndices(i),
        curve_descriptions(cd),
        totalCurveLength(tcl),
        smoothnessWeight(sw) {}
};

class Optimizer
{
public:
    Optimizer(int size);
    ~Optimizer();

    static void multiplyByA(const Vector& x, Vector &resultX, Vector &diagM, ProblemSettings &prob);

    static void multiplyByA(const Vector& x, const Vector& y, Vector &resultX, Vector &resultY, Grid &grid,
			    const std::vector<int> &curveIndices,
                            std::vector< std::vector<Intersection> > &mapCurveToConstraints,
                            float totalCurveLength, float smoothnessWeight);

    static void multiplyByAWithoutWeights(Vector& x, Vector& y, Vector &resultX,
                                          Vector &resultY, Grid &grid, std::vector<int> &curveIndices,
                                          std::vector< std::vector<Intersection> > &mapCurveToConstraints);

    void optimizeImplicitFastWithWeights(Grid &grid, int numberOfVectorFields,
                                         std::vector<PolygonalPath> curves,
                                         std::vector< std::pair<Vector*, Vector*> >& vectorFields,
                                         unsigned short *mapCurveToVectorField,
                                         float *mapCurveToError,
                                         float smoothnessWeight = 0.5);
};

#endif // OPTIMIZER_H
