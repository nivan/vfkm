#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <sstream>

#include "Optimizer.h"
#include "Util.h"

using namespace std;

void initExperiment(string filename, Cluster*& rootCluster, Grid*& g, 
		    vector<PolygonalPath>& curves, int gridResolution){
    //load curves
    float xmin, xmax, ymin, ymax, tmin, tmax;
    Util::loadCurves(filename, curves, xmin, xmax, ymin, ymax, tmin, tmax);

    //create grid
    g = new Grid(xmin,ymin,xmax-xmin,ymax-ymin,gridResolution, gridResolution);

    //create root cluster
    rootCluster = new Cluster;
    stringstream ss;
    ss << curves.size();
    rootCluster->name = ss.str();
    rootCluster->parent = NULL;
    Vector* rootVFX = new Vector(gridResolution * gridResolution);
    rootVFX->setValues(0.0);
    Vector* rootVFY = new Vector(gridResolution * gridResolution);
    rootVFY->setValues(0.0);

    rootCluster->vectorField = make_pair(rootVFX, rootVFY);

    for(size_t i = 0 ; i < curves.size() ; ++i){
        rootCluster->indices.push_back(i);
        rootCluster->curveErrors.push_back(0.0f);
    }
}

int main(int argc, char *argv[]){
    int rightNumberOfParameters = 6;

    if(argc != rightNumberOfParameters){
        //print usage
        cout << "./vfkm trajectoryFile gridResolution numberOfVectorFields smoothnessWeight" << endl;
        return 0;
    }

    //set parameters
    string filename(argv[1]);
    int    gridResolution = atoi(argv[2]);
    int    maxNumberOfVectorFields = atoi(argv[3]);
    float  smoothnessWeight = atof(argv[4]);

    //load files
    cout << "Loading Files..." << endl;
    vector<PolygonalPath> curves;
    Cluster* rootCluster = NULL;
    Grid* g = NULL;

    //for now always project
    cout << "Loading data" << endl;
    initExperiment(filename, rootCluster, g, curves, gridResolution);

    //optimize
    cout << "Optimizing..." << endl;
    float totalErrorsPerExecution[maxNumberOfVectorFields];
    for(int k = 1 ; k <= maxNumberOfVectorFields ; ++k){
	//
	Cluster* currentCluster = rootCluster;

	int numberOfVectorFields = k;


	//Optimize
	cout << "Optimizing " << k << " vector fields" << endl;
	Optimizer op(g->getResolutionX() * g->getResolutionY());
	int numberOfCurves = currentCluster->indices.size();
	unsigned short mapCurveToVF[numberOfCurves];
	float mapCurveToError[numberOfCurves];
	unsigned int mapCurveToIndexInCurveVector[numberOfCurves];

	vector<PolygonalPath> curvesInCurrentCluster;

	for(int i = 0 ; i <numberOfCurves ; ++i){
	    mapCurveToError[i] = 0;
	    mapCurveToVF[i] = -1;	   
	    curvesInCurrentCluster.push_back(curves.at(currentCluster->indices.at(i)));
	    mapCurveToIndexInCurveVector[i] = currentCluster->indices.at(i);
	}

	//create vector fields
	vector<pair<Vector*,Vector*> > vectorFields;
	int gridDimension  = g->getResolutionX() * g->getResolutionY();
	for(int i = 0 ; i < numberOfVectorFields ; ++i){
	    Vector* xComponent = new Vector(gridDimension);
	    Vector* yComponent = new Vector(gridDimension);
	    vectorFields.push_back(make_pair(xComponent, yComponent));
	}

	op.optimizeImplicitFastWithWeights(*g,numberOfVectorFields, curvesInCurrentCluster,
                                       vectorFields, &(mapCurveToVF[0]), mapCurveToError,smoothnessWeight);

	double totalError = 0.0f;

	for(int i = 0 ; i < numberOfCurves ; ++i){
	    totalError += mapCurveToError[i];	    
	}

	totalErrorsPerExecution[k-1] = totalError;

	//clean memory
	cout << "Cleaning Memory" << endl;
	for(int i = 0 ; i < numberOfVectorFields ; ++i){
	    std::pair<Vector*, Vector*> &vfs = vectorFields.at(i);
	    Vector* componentX  = vfs.first;
	    if(componentX != NULL){
		delete componentX;
	    }
	    
	    Vector* componentY  = vfs.second;
	    if(componentY != NULL){
		delete componentY;
	    }

	}
    }

    cout << endl << "Error Per Execution" << endl;
    for(int k = 1 ; k <= maxNumberOfVectorFields ; ++k){
	cout << "   With " << k << " vfs error was " << totalErrorsPerExecution[k-1] << endl;
    }

}
