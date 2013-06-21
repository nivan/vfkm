#include "Optimizer.h"
#include <PolygonalPath.h>
#include <Grid.h>

#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <cassert>

#include "ConstraintMatrix.h"

using namespace std;

#define xDEBUG
#define xWRITE_OUT_TO_FILE

Vector *Ax, *Ay;

/******************************************************************************/

Optimizer::Optimizer(int size)
{
  Ax	= new Vector(size);
  Ay	= new Vector(size);
}

Optimizer::~Optimizer()
{
  delete Ax;
  delete Ay;
}

double computeErrorImplicit
(const Grid &,
 const Vector &vfXComponent, const Vector& vfYComponent,
 float totalCurveLength,
 float smoothnessWeight,
 const CurveDescription &curve)
{
    double error = 0.0;
    Vector vx(2*curve.segments.size());
    Vector vy(2*curve.segments.size());
    for (size_t i = 0; i<curve.segments.size(); ++i) {
        curve.segments[i].add_cx(vx, vfXComponent);
        curve.segments[i].add_cx(vy, vfYComponent);
    }

    vx -= curve.rhsx;
    vy -= curve.rhsy;

    // LT . L = [[1/3 1/6] [1/6 1/3]]
    for (int i=0; i<vx.getDimension(); i+=2) {
        double this_error_x = (vx[i] * vx[i] + vx[i] * vx[i+1] + vx[i+1] * vx[i+1]) / 3.0;
        double this_error_y = (vy[i] * vy[i] + vy[i] * vy[i+1] + vy[i+1] * vy[i+1]) / 3.0;
        error += (this_error_x + this_error_y) * curve.length;
    }

    assert(error >= 0.0);
    return error * (1.0 - smoothnessWeight) / (totalCurveLength);
}

//void computeErrorImplicitConstraints(Grid &g, const Vector &vfXComponent, const Vector& vfYComponent,
//                                     int curveIndex, vector<vector<Intersection> >& mapCurveToConstraints){

//    vector<Intersection> &constraints = mapCurveToConstraints[curveIndex];
//    int numberOfConstraints = constraints.size();

//    Vector cx(3 * (numberOfConstraints / 2));
//    cx.setValues(0.0f);

//    for(int i = 0 ; i < numberOfConstraints ; i += 2){
//       Intersection currentIntersection = constraints.at(i);
//       Intersection nextIntersection = constraints.at(i+1);

//       //cx[vfXComponent[currentIntersection.indexV1];
//    }
//}

//////// With Weights
#define xDEBUG
#define EPSILON 0.000000001

template <int P>
struct FastPow
{
    static inline float value(float x) { return x * FastPow<P-1>::value(x); }
};

template <>
struct FastPow<1>
{
    static inline float value(float x) { return x; }
};

void Optimizer::multiplyByA(const Vector& x, Vector &resultX, Vector &diagATA, ProblemSettings &prob)
{
    Grid &grid = prob.grid;
    const vector<int> &curveIndices = prob.curveIndices;
    const vector<CurveDescription> &curve_descriptions = prob.curve_descriptions;
    VECTOR_TYPE totalCurveLength = prob.totalCurveLength;
    VECTOR_TYPE smoothnessWeight = prob.smoothnessWeight;

    Vector Ax(x);
    resultX.setValues(0.0);

    int numberOfVertices = grid.getResolutionX() * grid.getResolutionY();

    // L . x
    grid.multiplyByLaplacian2(Ax, diagATA); // gets overwritten the second time, but whatever.
    // L^T . L . x
    grid.multiplyByLaplacian2(Ax, diagATA);
    Ax.scale(smoothnessWeight/numberOfVertices);

    for (size_t k=0; k<curveIndices.size(); ++k) {
        int i = curveIndices[k];
        const CurveDescription &curve = curve_descriptions[i];
        float k_global = (1.0f - smoothnessWeight) / totalCurveLength;
        curve.add_cTcx(resultX, x, k_global);
    }

    resultX.add(Ax);
}

/******************************************************************************/

void cg_solve(ProblemSettings &prob, const Vector &b, Vector &x)
{
    // mat m(prob);
    // DiagonalPrec prec(m);
    int max_iter = 10000;
    VECTOR_TYPE tol = 1e-8;
    VECTOR_TYPE resid;
    // CG(m, x, b, prec, max_iter, tol);

    VECTOR_TYPE normb = b.length();
    Vector r(x.getDimension());
    Vector z(x.getDimension());
    Vector q(x.getDimension());
    Vector p(x.getDimension());
    Vector diagATA(x.getDimension());

    VECTOR_TYPE alpha, beta, rho = 1, rho_1 = 1;

    diagATA.setValues(0.0);
    Vector t(x.getDimension());
    Optimizer::multiplyByA(x, t, diagATA, prob);

    r.setValues(b);
    r -= t;
    if (normb == 0.0) {
        normb = 1.0;
    }    
    resid = r.length() / normb;

    if (resid <= tol) {
        tol = resid;
        max_iter = 0;
    }
    for (int i=1; i<=max_iter; ++i) {
        z.setValues(r);
        // z/=diagATA; // precondition

        rho = r.dot(z);

        if (i == 1) {
            p.setValues(z);
        } else {
            beta = rho / rho_1;
            p.add_scale(z, 1.0, beta);
        }

        Optimizer::multiplyByA(p, q, t, prob);
        alpha = rho / p.dot(q);
        x.add_scale(p,  alpha);
        r.add_scale(q, -alpha);

        resid = r.length() / normb;

        //cout << "Resid " << resid << endl;

        if (resid <= tol) {
            tol = resid;
            max_iter = i;
            break;
        }
        rho_1 = rho;
    }
};

void print_problem(ProblemSettings &prob, const Vector &bx, const Vector &by)
{
    int v = prob.grid.getResolutionX() * prob.grid.getResolutionY();
    Vector foo(v);
    for (int i=0; i<v; ++i) {
        Vector t(v), output(v);
        t[i] = 1;
        Optimizer::multiplyByA(t, output, foo, prob);
        cout << output.toString() << endl;
    }
    cout << bx.toString() << endl;
    cout << by.toString() << endl;
}

void optimizeVectorFieldWithWeights(Grid &grid, Vector &initialGuessX, Vector &initialGuessY,
                                    const vector<int> &curveIndices,
                                    const vector<CurveDescription> &curve_descriptions,
                                    float totalCurveLength,
                                    float smoothnessWeight)
{
    // cout << "Fitting Vector Field" << endl;
    ///compute independent terms
    int numberOfVertices = grid.getResolutionX() * grid.getResolutionY();

    Vector 
        indepx(numberOfVertices), 
        indepy(numberOfVertices);
    indepx.setValues(0.0f);
    indepy.setValues(0.0f);

    for(size_t k = 0; k < curveIndices.size() ; ++k) {
        int i = curveIndices[k];
        const CurveDescription &curve = curve_descriptions[i];

        for (size_t j=0; j<curve.segments.size(); ++j) {
            float k = (1.0 - smoothnessWeight) * (curve.segments[j].time[1] - curve.segments[j].time[0])/totalCurveLength;
            curve.segments[j].add_cTx(indepx, curve.rhsx, k);
            curve.segments[j].add_cTx(indepy, curve.rhsy, k);
        }
    }

    ProblemSettings prob(grid, curveIndices, curve_descriptions,
                         totalCurveLength, smoothnessWeight);

    Vector x(initialGuessX), y(initialGuessY);

    cg_solve(prob, indepx, x);
    cg_solve(prob, indepy, y);
    initialGuessX.setValues(x);
    initialGuessY.setValues(y);
}

pair< vector<int>, vector< vector<int> > > compute_first_assignment
(Grid &grid, int numberOfVectorFields,
 const vector<CurveDescription> &curves,
 float totalCurveLength,
 float smoothnessWeight)
{
    vector<double> errors(curves.size(), 1e10);

    vector<pair<Vector, Vector> > vector_fields(numberOfVectorFields);

    for (int i=0; i<numberOfVectorFields; ++i) {
        vector<int> curveIndices;
//        if(i == 0){
//            srand(time(NULL));
//            int index = ((1.0f*rand())/RAND_MAX) * curves.size();
//            cout << "Initial index " << index << endl;
//            curveIndices.push_back(index);
//        }
//        else
        curveIndices.push_back(max_element(errors.begin(), errors.end()) - errors.begin());

        //optimize
        pair<Vector, Vector> &vs = vector_fields[i];
        vs.first = Vector(grid.getResolutionX() * grid.getResolutionY());
        vs.second = Vector(grid.getResolutionX() * grid.getResolutionY());
        Vector &xComponent(vs.first);
        Vector &yComponent(vs.second);
        xComponent.setValues(0.0);
        yComponent.setValues(0.0);

        optimizeVectorFieldWithWeights(grid, xComponent, yComponent,
                                       curveIndices, curves, totalCurveLength, smoothnessWeight);

        for (size_t j=0; j < curves.size(); ++j) {
            errors[j] = min(errors[j], computeErrorImplicit(grid, xComponent, yComponent, totalCurveLength, smoothnessWeight, curves[j]));
        }
    }

    vector<int> result;
    vector<vector<int> > result_indices(numberOfVectorFields);

    for (size_t i=0; i<curves.size(); ++i) {
        vector<double> curve_errors;
        for (int j=0; j<numberOfVectorFields; ++j)
            curve_errors.push_back(computeErrorImplicit(grid, vector_fields[j].first, vector_fields[j].second, totalCurveLength, smoothnessWeight, curves[i]));
        int best_index = min_element(curve_errors.begin(), curve_errors.end()) - curve_errors.begin();
        result_indices[best_index].push_back(result.size());
        result.push_back(best_index);
    }
    return make_pair(result, result_indices);
}

void set_constraints(vector<CurveDescription> &curve_descriptions,
                     float &totalCurveLength,
                     vector<PolygonalPath> &curves,                     
                     const Grid &grid)
{
    totalCurveLength = 0.0f;
    int numberOfCurves = curves.size();


    for(int i = 0 ; i < numberOfCurves ; ++i) {

        bool bad_break = false;        

        PolygonalPath &p = curves.at(i);

        for (size_t j=0; j<p.numberOfPoints()-1; ++j) {
            if (p.getPoint(j+1).second < p.getPoint(j).second) {
                cerr << "Line is broken, has backward time." << endl;
//                exit(1);
            }
        }

        grid.clipLine(p);
        for (size_t j=0; j<p.numberOfPoints()-1; ++j) {
            if (p.getPoint(j+1).second < p.getPoint(j).second) {
                cerr << i << " - Line clipper is broken, introduced backward time: ";
                cerr << p.getPoint(j+1).second << " " << p.getPoint(j).second << endl;
                // exit(1);
                bad_break = true;
                break;
            }
        }

        CurveDescription curve;
        curve.index = i;
        curve.length = 0;
        if (!bad_break) {
            curve = grid.curve_description(p);
            totalCurveLength += curve.length;
        }
        curve_descriptions.push_back(curve);
    }
}

typedef vector< pair<Vector, Vector> > AllVectorFields;
typedef vector< CurveDescription > AllConstraints;

void optimize_all_vector_fields(
        AllVectorFields &vectorFields,
        Grid &grid,
        const vector<vector<int> > &mapVectorFieldCurves,
        const AllConstraints &curves,
        float totalCurveLength,
        float smoothnessWeight)
{
    size_t numberOfVectorFields = vectorFields.size();
    //optimize each vector field
    for(size_t j = 0 ; j < numberOfVectorFields ; ++j) {
        pair<Vector, Vector> &currentVectorField = vectorFields.at(j);
        const vector<int>& curveIndices = mapVectorFieldCurves.at(j);

        //optimize
        Vector &xComponent = currentVectorField.first;
        Vector &yComponent = currentVectorField.second;

        optimizeVectorFieldWithWeights
            (grid, xComponent, yComponent,
             curveIndices, curves,
             totalCurveLength, smoothnessWeight);
    }
}

double get_total_error(const vector<CurveDescription> &curves,
                       const AllVectorFields &vectorFields,
                       const unsigned short *mapCurveToVectorField,
                       float totalCurveLength,
                       float smoothnessWeight,
                       const Grid &grid)
{
    double totalError = 0.0;
    size_t numberOfCurves = curves.size();
    int numberOfVectorFields = vectorFields.size();
    vector<float> lengths(numberOfVectorFields, 0.0f);
    for (size_t i = 0; i < numberOfCurves; ++i) {
        const CurveDescription& currentCurve = curves.at(i);
        int vectorFieldIndex = mapCurveToVectorField[i];
        const pair<Vector, Vector> &currentVectorField = vectorFields.at(vectorFieldIndex);
        double error = computeErrorImplicit(grid, currentVectorField.first, currentVectorField.second, totalCurveLength, smoothnessWeight, currentCurve);
        totalError += error;
        lengths[vectorFieldIndex] += currentCurve.length;
    }

    for (int i=0; i < numberOfVectorFields; ++i) {
        const pair<Vector, Vector> &currentVectorField = vectorFields.at(i);
        Vector t1(currentVectorField.first), t2(currentVectorField.second);
        grid.multiplyByLaplacian(t1, t2);
        totalError += t1.length2() * smoothnessWeight * (lengths[i] / totalCurveLength);
        totalError += t2.length2() * smoothnessWeight * (lengths[i] / totalCurveLength);
    }

    return totalError;
}

void optimize_assignments(int &total_change,
                          double &totalError,
                          unsigned short *mapCurveToVectorField,
                          vector<vector<int> > &mapVectorFieldCurves,
                          float *mapCurveToError,
                          const AllVectorFields &vectorFields,
                          const vector<CurveDescription> &curves,
                          float totalCurveLength,
                          float smoothnessWeight,
                          Grid &grid)
{
    //updating mapVectorFieldCurves
    totalError = 0.0;
    total_change = 0;
    size_t numberOfCurves = curves.size();
    size_t numberOfVectorFields = vectorFields.size();
    for(size_t i = 0 ; i < numberOfCurves ; ++i){
        bool change = false;
        const CurveDescription& currentCurve = curves.at(i);

        size_t vectorFieldIndex = mapCurveToVectorField[i];
        int newVectorFieldIndex = vectorFieldIndex;
        const pair<Vector, Vector> &currentVectorField = vectorFields.at(vectorFieldIndex);
        double error = computeErrorImplicit(grid, currentVectorField.first, currentVectorField.second, totalCurveLength, smoothnessWeight, currentCurve);

        for(size_t j = 0 ; j < numberOfVectorFields ; ++j){
            if(j == vectorFieldIndex)
                continue;

            const pair<Vector, Vector> &vectorField = vectorFields.at(j);
            double currentError = computeErrorImplicit(grid, vectorField.first, vectorField.second, totalCurveLength, smoothnessWeight, currentCurve);

            if(currentError < error){
                newVectorFieldIndex = j;
                error = currentError;
                change = true;
            }
        }
        total_change += change;

        totalError += error;

        //assign
        mapCurveToVectorField[i] = newVectorFieldIndex;
        mapCurveToError[i] = error;
    }

    //updating mapVectorFieldCurves
    for(size_t i = 0 ; i < numberOfVectorFields ; ++i){
        vector<int> &container = mapVectorFieldCurves.at(i);
        container.clear();
    }

    for(size_t i = 0 ; i < numberOfCurves ; ++i){
        int vectorFieldIndex = mapCurveToVectorField[i];
        vector<int> &vectorFieldCurves = mapVectorFieldCurves.at(vectorFieldIndex);
        vectorFieldCurves.push_back(i);
    }
}

void repopulate_empty_cluster(vector<vector<int> > &mapVectorFieldCurves,
                              unsigned short *mapCurveToVectorField,
                              AllVectorFields &vectorFields)
{
    size_t numberOfVectorFields = vectorFields.size();

    for (size_t i=0; i<numberOfVectorFields; ++i) {
        vector<int> &container = mapVectorFieldCurves.at(i);
        if (container.size() == 0) {
            vectorFields[i].first.setValues(0.0f);
            vectorFields[i].second.setValues(0.0f);
            int max_index = -1;
            size_t sz = 0;
            for (size_t j=0; j<numberOfVectorFields; ++j) {
                if (mapVectorFieldCurves[j].size() > sz) {
                    sz = mapVectorFieldCurves[j].size();
                    max_index = j;
                }
            }
            vector<int> n1, n2;
            for (size_t j=0; j<sz; ++j) {
                int curve = mapVectorFieldCurves[max_index][j];
                if (j % 2) {
                    n1.push_back(curve);
                    mapCurveToVectorField[curve] = i;
                } else {
                    n2.push_back(curve);
                    mapCurveToVectorField[curve] = max_index;
                }
            }
            mapVectorFieldCurves[i] = n1;
            mapVectorFieldCurves[max_index] = n2;
        }
    }
}

void Optimizer::optimizeImplicitFastWithWeights
    (Grid &grid, int numberOfVectorFields,
     vector<PolygonalPath> curves,
     std::vector< pair<Vector*, Vector*> >& finalVectorFields,
     unsigned short *mapCurveToVectorField,
     float *mapCurveToError,
     float smoothnessWeight)
{
    float totalCurveLength;
    vector<CurveDescription> curve_descriptions;

    // create vector fields
    int sz = grid.getResolutionX() * grid.getResolutionY();
    vector< pair<Vector, Vector> > vectorFields(numberOfVectorFields, make_pair(Vector(sz), Vector(sz)));

    // determine constraints
    set_constraints(curve_descriptions, totalCurveLength, curves, grid);
    // first assignment
    pair<vector<int>, vector<vector<int> > > f = compute_first_assignment(
        grid, numberOfVectorFields, curve_descriptions, totalCurveLength, smoothnessWeight);
    vector<vector<int> > mapVectorFieldCurves = f.second;
    copy(f.first.begin(), f.first.end(), mapCurveToVectorField);

    //optimize
    int numberOfIterations = 0;
    double totalError = 1e20;

    while(numberOfIterations < 100){
        int total_change = 0;

        cout << "Before optimization: " << totalError << endl;
        optimize_all_vector_fields(vectorFields, grid, mapVectorFieldCurves, curve_descriptions, totalCurveLength, smoothnessWeight);

        totalError = get_total_error(curve_descriptions, vectorFields, mapCurveToVectorField, totalCurveLength, smoothnessWeight, grid);
        cout << "After optimization: " << totalError << endl;

        optimize_assignments(total_change, totalError, mapCurveToVectorField, mapVectorFieldCurves, mapCurveToError, vectorFields, curve_descriptions, totalCurveLength, smoothnessWeight, grid);
        totalError = get_total_error(curve_descriptions, vectorFields, mapCurveToVectorField, totalCurveLength, smoothnessWeight, grid);

        cout << "After assignment: " << totalError << " changes: " << total_change << endl;

        repopulate_empty_cluster(mapVectorFieldCurves, mapCurveToVectorField, vectorFields);

        ++numberOfIterations;
        if(total_change == 0)
            break;
    }

    for(int i = 0 ; i < numberOfVectorFields ; ++i){
        finalVectorFields[i].first->setValues(vectorFields[i].first);
        finalVectorFields[i].second->setValues(vectorFields[i].second);
    }
}

