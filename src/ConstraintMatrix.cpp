#include "ConstraintMatrix.h"
#include <cassert>
#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>

using namespace std;

#define xDEBUG

ConstraintMatrix::ConstraintMatrix()
{
}

void ConstraintMatrix::multiply(const std::vector<Intersection>& constraints, Vector &v){

    cout << "TO DO!!!:multiply" << endl;
    exit(1);

    Vector result(constraints.size());

    vector<Intersection>::const_iterator iterator;

    result.setValues(0.0);

    for(iterator = constraints.begin() ; iterator != constraints.end() ; ++iterator){
        Intersection inter = *iterator;

        float v1 = v[inter.indexV1];
        v1 *= (1.0 - inter.lambda);

        result.setValue(inter.indexV1, v1);

        float v2 = v[inter.indexV2];
        v2 *= inter.lambda;

        result.setValue(inter.indexV2, v2);
    }

    v.setValues(result);
}

void ConstraintMatrix::multiply(const vector<Intersection>& constraints, const Vector &v, Vector &result){

    assert(constraints.size() == (size_t) result.getDimension());

    vector<Intersection>::const_iterator iterator;

    result.setValues(0.0);

    int countConstraint = 0;

    for(iterator = constraints.begin() ; iterator != constraints.end() ; ++iterator){
        Intersection inter = *iterator;

        float v1 = v[inter.indexV1];
        v1 *= (1.0 - inter.lambda);

        float v2 = v[inter.indexV2];
        v2 *= inter.lambda;

        result.setValue(countConstraint, v1 + v2);

        ++countConstraint;
    }

}

void ConstraintMatrix::multiply(const std::vector<Intersection>& constraints, const Vector &firstComponent, const Vector &secondComponent,
                                Vector &firstResult, Vector &secondResult){
    assert(constraints.size() == (size_t) firstResult.getDimension() &&
           constraints.size() == (size_t) secondResult.getDimension());

    vector<Intersection>::const_iterator iterator;

    //result.setValues(0.0);
    firstResult.setValues(0.0);
    secondResult.setValues(0.0);

    int countConstraint = 0;

    for(iterator = constraints.begin() ; iterator != constraints.end() ; ++iterator){
        Intersection inter = *iterator;

        //first
        float firstv1 = firstComponent[inter.indexV1];
        firstv1 *= (1.0 - inter.lambda);

        float firstv2 = firstComponent[inter.indexV2];
        firstv2 *= inter.lambda;

        firstResult.setValue(countConstraint, firstv1 + firstv2);

        //second
        float secondv1 = secondComponent[inter.indexV1];
        secondv1 *= (1.0 - inter.lambda);

        float secondv2 = secondComponent[inter.indexV2];
        secondv2 *= inter.lambda;

        secondResult.setValue(countConstraint, secondv1 + secondv2);

        ++countConstraint;
    }
}

void ConstraintMatrix:: multiplyTranspose(const std::vector<Intersection>& constraints, std::map<int,vector<pair<int,Intersection> > >& mapVertexConstraint,
                                          const Vector &v, Vector &result){
    assert(constraints.size() == (size_t) v.getDimension());

    int numberOfVertices = result.getDimension();

        //multiply
        for(int i = 0 ; i < numberOfVertices ; ++i){
            vector<pair<int,Intersection> > &v1 = mapVertexConstraint[i];

            int numberOfConstraintsForV1 = v1.size();


#ifdef DEBUG
            cout << "Number of constrains for vertex " << i << " is " << numberOfConstraintsForV1 << endl;
#endif


            float finalValue = 0;

            for(int j = 0 ; j < numberOfConstraintsForV1 ; ++j){
                pair<int,Intersection> indexInter = v1.at(j);
                Intersection inter = indexInter.second;

#ifdef DEBUG
                cout << "Intersection " << j << " is" << endl;
                cout << "Index v1 " << inter.indexV1 << " index v2 " << inter.indexV2 << " "
                     << "lambda " << inter.lambda << " direction " << inter.dirAtIntersection.toString() << endl;
#endif

                //decide the right coefficient
                float coef;

                if(inter.indexV1 == i){
                    coef = 1.0 - inter.lambda;
                }
                else if(inter.indexV2 == i){
                    coef = inter.lambda;
                }
                else{
                    cout << "Error while multiplying: Invalid mapping" << endl;
                    exit(1);
                }

                //multiply
                finalValue += (coef * v[indexInter.first]);
            }

            //cout << "Final Value " << finalValue << endl;
            result.setValue(i, finalValue);
        }

}

void ConstraintMatrix::multiplyTranspose(vector<Intersection>& constraints, const Vector &v, Vector &result){    
    int numberOfConstraints = constraints.size();

    assert(numberOfConstraints == v.getDimension());

    //multiply
    result.setValues(0.0);
    for(int i = 0 ; i < numberOfConstraints ; ++i){
        Intersection inter = constraints.at(i);
        float vi = v[i];

        result.setValue(inter.indexV1, (vi * (1.0 - inter.lambda)) + result[inter.indexV1]);
        result.setValue(inter.indexV2, (vi * inter.lambda) + result[inter.indexV2]);
    }

}

void ConstraintMatrix::multiplyTranspose(vector<Intersection>& constraints, const Vector &v,
                                         const Vector &w, Vector &resultV, Vector &resultW){
    int numberOfVertices = resultV.getDimension();
    int numberOfConstraints = constraints.size();

    assert(numberOfConstraints == v.getDimension()&&
           numberOfConstraints == w.getDimension() && resultW.getDimension() == numberOfVertices);

    //multiply
    resultV.setValues(0.0);
    resultW.setValues(0.0);
    for(int i = 0 ; i < numberOfConstraints ; ++i){
        Intersection inter = constraints.at(i);

        //V
        float vi = v[i];
        resultV.setValue(inter.indexV1, (vi * (1.0 - inter.lambda)) + resultV[inter.indexV1]);
        resultV.setValue(inter.indexV2, (vi * inter.lambda) + resultV[inter.indexV2]);

        //W
        float wi = w[i];
        resultW.setValue(inter.indexV1, (wi * (1.0 - inter.lambda)) + resultW[inter.indexV1]);
        resultW.setValue(inter.indexV2, (wi * inter.lambda) + resultW[inter.indexV2]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

//void ConstraintMatrix:: multiplyTranspose(const std::vector<Intersection>& constraints, std::map<int,vector<pair<int,Intersection> > >& mapVertexConstraint,
//                                         const Vector &firstComponent, const Vector &secondComponent,
//                                         Vector &firstResult, Vector &secondResult){
//    assert(constraints.size() == (size_t) firstComponent.getDimension() &&
//           constraints.size() == (size_t) secondComponent.getDimension() &&
//           firstResult.getDimension() == secondResult.getDimension());

//    int numberOfVertices = firstResult.getDimension();

//    //multiply
//    for(int i = 0 ; i < numberOfVertices ; ++i){
//        vector<pair<int,Intersection> > &v1 = mapVertexConstraint[i];

//        int numberOfConstraintsForV1 = v1.size();

//#ifdef DEBUG
//        cout << "Number of constrains for vertex " << i << " is " << numberOfConstraintsForV1 << endl;
//#endif


//        float finalValue1 = 0;
//        float finalValue2 = 0;

//        for(int j = 0 ; j < numberOfConstraintsForV1 ; ++j){
//            pair<int,Intersection> indexInter = v1.at(j);
//            Intersection inter = indexInter.second;

//#ifdef DEBUG
//            cout << "Intersection " << j << " is" << endl;
//            cout << "Index v1 " << inter.indexV1 << " index v2 " << inter.indexV2 << " "
//                 << "lambda " << inter.lambda << " direction " << inter.dirAtIntersection.toString() << endl;
// #endif

//            //decide the right coefficient
//            float coef;

//            if(inter.indexV1 == i){
//                coef = 1.0 - inter.lambda;
//            }
//            else if(inter.indexV2 == i){
//                coef = inter.lambda;
//            }
//            else{
//                cout << "Error while multiplying: Invalid mapping" << endl;
//                exit(1);
//            }

//            //multiply
//            finalValue1 += (coef * firstComponent[indexInter.first]);
//            finalValue2 += (coef * secondComponent[indexInter.first]);
//        }

//        //cout << "Final Value " << finalValue << endl;
//        firstResult.setValue(i, finalValue1);
//        secondResult.setValue(i, finalValue2);
//    }
//}

//void ConstraintMatrix::multiplyTranspose(vector<Intersection>& constraints, const Vector &v, Vector &result){
//    assert(constraints.size() == v.getDimension());

//    int numberOfConstraints = constraints.size();
//    int numberOfVertices = result.getDimension();

//    map<int, vector<pair<int,Intersection> > > mapVertexConstraint;

//    //build map vertex constraint
//    for(int i = 0 ; i < numberOfVertices ; ++i){
//        mapVertexConstraint[i] = vector<pair<int,Intersection> >();
//    }

//    for(int i = 0 ; i < numberOfConstraints ; ++i){
//        Intersection inter = constraints.at(i);
//        vector<pair<int, Intersection> > &v1 = mapVertexConstraint[inter.indexV1];
//        vector<pair<int, Intersection> > &v2 = mapVertexConstraint[inter.indexV2];

//        v1.push_back(make_pair(i, inter));
//        v2.push_back(make_pair(i, inter));
//    }

//    //multiply
//    for(int i = 0 ; i < numberOfVertices ; ++i){
//        vector<pair<int,Intersection> > &v1 = mapVertexConstraint[i];

//        int numberOfConstraintsForV1 = v1.size();

//        float finalValue = 0;

//        for(int j = 0 ; j < numberOfConstraintsForV1 ; ++j){
//            pair<int,Intersection> indexInter = v1.at(j);
//            Intersection inter = indexInter.second;

//            //decide the right coefficient
//            float coef;

//            if(inter.indexV1 == i){
//                coef = 1.0 - inter.lambda;
//            }
//            else if(inter.indexV2 == i){
//                coef = inter.lambda;
//            }
//            else{
//                cout << "Error while multiplying: Invalid mapping" << endl;
//                exit(1);
//            }

//            //multiply
//            finalValue += (coef * v[indexInter.first]);
//        }

//        //cout << "Final Value " << finalValue << endl;
//        result.setValue(i, finalValue);
//    }
//}

///////////////////////////////////////////////////////////////////////////////////////////////

void ConstraintMatrix::multiplyCTC(vector<Intersection>& constraints, const Vector &v, Vector &result){
    assert(v.getDimension() == result.getDimension());

    Vector result1(constraints.size());
    ConstraintMatrix::multiply(constraints, v, result1);
    ConstraintMatrix::multiplyTranspose(constraints, result1, result);
}

void ConstraintMatrix::multiplyCTC(std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& mapVertexConstraint, const Vector &v, Vector &result){
    assert(v.getDimension() == result.getDimension());

    Vector result1(constraints.size());

    ConstraintMatrix::multiply(constraints, v, result1);
    ConstraintMatrix::multiplyTranspose(constraints, mapVertexConstraint, result1, result);
}

void ConstraintMatrix::multiplyCTC(std::vector<Intersection>& constraints, std::map<int,std::vector<std::pair<int,Intersection> > >& ,
                                   const Vector &firstComponent, const Vector &secondComponent,
                                   Vector &, Vector &secondResult){
    exit(1);
    assert(firstComponent.getDimension() == firstComponent.getDimension() &&
           secondComponent.getDimension() == secondResult.getDimension());

    Vector result1(constraints.size());
    Vector result2(constraints.size());

    ConstraintMatrix::multiply(constraints, firstComponent, secondComponent, result1, result2);
    //ConstraintMatrix::multiplyTranspose(constraints, mapVertexConstraint, result1, result2, firstResult, secondResult);
}

void ConstraintMatrix::multiplyCTC(std::vector<Intersection>& constraints,
                                   const Vector &firstComponent, const Vector &secondComponent,
                                   Vector &firstResult, Vector &secondResult){
    assert(firstComponent.getDimension() == firstComponent.getDimension() &&
           secondComponent.getDimension() == secondResult.getDimension());

    Vector result1(constraints.size());
    Vector result2(constraints.size());

    ConstraintMatrix::multiply(constraints, firstComponent, secondComponent, result1, result2);
    ConstraintMatrix::multiplyTranspose(constraints, result1, result2, firstResult, secondResult);
}
