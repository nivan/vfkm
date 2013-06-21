#ifndef CONSTRAINTMATRIX_H
#define CONSTRAINTMATRIX_H
#include <map>
#include "Grid.h"



class ConstraintMatrix
{
public:
    ConstraintMatrix();

    static void multiply(const std::vector<Intersection>& constraints, const Vector &, Vector &result);
    static void multiply(const std::vector<Intersection>& constraints, Vector &);
    static void multiply(const std::vector<Intersection>& constraints, const Vector &, const Vector&, Vector &, Vector&);

    ///////////
    static void multiplyTranspose(const std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& mapVertexConstraint,
                                  const Vector &, Vector &result);
    static void multiplyTranspose(std::vector<Intersection>& constraints, const Vector &, Vector &result);

    static void multiplyTranspose(std::vector<Intersection>& constraints, const Vector &v, const Vector &w, Vector &resultV, Vector &resultW);

    static void multiplyTranspose(std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& mapVertexConstraint,
                                  const Vector &firstComponent, const Vector &secondComponent,
                                  Vector &firstResult, Vector &secondResult, float weight);   

//    static void multiplyTranspose(std::vector<Intersection>& constraints, const Vector &v,
//                                  const Vector &w, Vector &resultV, Vector &resultW,
//                                  float weight);

    ///////////
    static void multiplyCTC(std::vector<Intersection>& constraints, const Vector &, Vector &result);
    static void multiplyCTC(std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& mapVertexConstraint, const Vector &, Vector &result);

    static void multiplyCTC(std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& mapVertexConstraint,
                            const Vector &firstComponent, const Vector &secondComponent,
                            Vector &firstResult, Vector &secondResult);

    static void multiplyCTC(std::vector<Intersection>& constraints,
                            const Vector &firstComponent, const Vector &secondComponent,
                            Vector &firstResult, Vector &secondResult);
};

#endif // CONSTRAINTMATRIX_H
