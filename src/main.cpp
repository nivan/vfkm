#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <map>
#include <cassert>

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

void saveExperiment(string directory, string currentFileLoaded, Cluster* root){
    stringstream ss;
    ss << directory << "/experiment.txt";

    ofstream experimentFile(ss.str().c_str());
    experimentFile << currentFileLoaded << endl;

    //write hierarchy
    queue<Cluster*> nodesToProcess;

    nodesToProcess.push(root);

    //hack
    map<Cluster*, string> mapClusterPath;
    mapClusterPath[root] = "r";

    experimentFile << "-1 r " << endl;

    while(!nodesToProcess.empty()){
        Cluster* c = nodesToProcess.front();
        nodesToProcess.pop();

        string name = mapClusterPath[c];

        //write curve indices file
        stringstream curveFileName;
        curveFileName << directory << "/curves_" << name << ".txt";
        ofstream curveIndicesFile(curveFileName.str().c_str());
        int numberOfCurves = c->indices.size();

        //cout << "NumCurves " << numberOfCurves << " errors " << c->curveErrors.size() << endl;
        assert(numberOfCurves == (int)c->curveErrors.size());

        for(int i = 0 ; i < numberOfCurves ; ++i){
            curveIndicesFile << c->indices.at(i) << " " << c->curveErrors.at(i) << endl;
        }
        curveIndicesFile.close();

//write vector field file file
        stringstream vectorFieldFileName;
        vectorFieldFileName << directory << "/vf_" << name << ".txt";
        ofstream vectorFieldFile(vectorFieldFileName.str().c_str());
        Vector* xComponent = c->vectorField.first;
        Vector* yComponent = c->vectorField.second;
        int gridDimension = xComponent->getDimension();

        vectorFieldFile << gridDimension << endl;

        for(int i = 0 ; i < gridDimension ; ++i){
            vectorFieldFile << xComponent[0][i] << " " << yComponent[0][i] << endl;
        }
        vectorFieldFile.close();

        //process children
        int numberOfChildren = c->children.size();

        for(int i = 0 ; i < numberOfChildren ; ++i){
            Cluster* child = c->children.at(i);

            stringstream ss;
            ss << name << "_" <<  i;
            mapClusterPath[child] = ss.str();

            experimentFile << name << " " << ss.str() << endl;

            nodesToProcess.push(child);
        }
    }

}

int main(int argc, char *argv[]){
    int rightNumberOfParameters = 6;

    if(argc != rightNumberOfParameters){
        //print usage
        cout << "./vfkm trajectoryFile gridResolution numberOfVectorFields smoothnessWeight outputDirectory" << endl;
        return 0;
    }

    //set parameters
    string filename(argv[1]);
    int    gridResolution = atoi(argv[2]);
    int    numberOfVectorFields = atoi(argv[3]);
    float  smoothnessWeight = atof(argv[4]);
    string outputDirectory(argv[5]);

    //load files
    cout << "Loading Files..." << endl;
    vector<PolygonalPath> curves;
    Cluster* rootCluster = NULL;
    Grid* g = NULL;

    //
    cout << "Loading data" << endl;
    initExperiment(filename, rootCluster, g, curves, gridResolution);

    //optimize
    cout << "Optimizing..." << endl;
    Cluster* currentCluster = rootCluster;

    //Optimize
    Optimizer op(g->getResolutionX() * g->getResolutionY());
    int numberOfCurves = currentCluster->indices.size();

    unsigned short mapCurveToVF[numberOfCurves];
    float mapCurveToError[numberOfCurves];
    unsigned int mapCurveToIndexInCurveVector[numberOfCurves];
    vector<float> mapVectorFieldToError;
    vector<PolygonalPath> curvesInCurrentCluster;

    for(int i = 0 ; i <numberOfCurves ; ++i){
        mapCurveToError[i] = 0;
        mapCurveToVF[i] = -1;

        curvesInCurrentCluster.push_back(curves.at(currentCluster->indices.at(i)));
        mapCurveToIndexInCurveVector[i] = currentCluster->indices.at(i);
    }

    vector<pair<Vector*,Vector*> > vectorFields;
    int gridDimension  = g->getResolutionX() * g->getResolutionY();
    for(int i = 0 ; i < numberOfVectorFields ; ++i){
        Vector* xComponent = new Vector(gridDimension);
        Vector* yComponent = new Vector(gridDimension);
        vectorFields.push_back(make_pair(xComponent, yComponent));
    }

    op.optimizeImplicitFastWithWeights(*g,numberOfVectorFields, curvesInCurrentCluster,
                                       vectorFields, &(mapCurveToVF[0]), mapCurveToError,smoothnessWeight);

    //count number of curves for each vf
    mapVectorFieldToError = vector<float>(vectorFields.size(),0);
    int numberOfChildren = vectorFields.size();
    int numberOfCurvesPerVF[numberOfChildren];
    for(int i = 0 ; i < numberOfChildren ; ++i){
        numberOfCurvesPerVF[i] = 0;
    }

    //compute error by vector field and total error
    vector<vector<int> > mapCurvesToClusters(numberOfChildren,vector<int>());
    vector<vector<float> > mapCurveErrorsToClusters(numberOfChildren,vector<float>());


    for(int i = 0 ; i < numberOfCurves ; ++i){
        vector<int>& curveCluster
                = mapCurvesToClusters.at(mapCurveToVF[i]);
        vector<float>& curveErrorsCluster
                = mapCurveErrorsToClusters.at(mapCurveToVF[i]);

        curveCluster.push_back(mapCurveToIndexInCurveVector[i]);
        curveErrorsCluster.push_back(mapCurveToError[i]);

        float &error = mapVectorFieldToError.at(mapCurveToVF[i]);
        error += mapCurveToError[i];
        numberOfCurvesPerVF[mapCurveToVF[i]] += 1;
    }

    //update cluster struct
    currentCluster->clearChildren();
    for(int i = 0 ; i < numberOfChildren ; ++i){
        Cluster* c = new Cluster();
        c->children.clear();
        c->parent = currentCluster;
        stringstream ss;
        ss << currentCluster->name << ":" << i;
        c->name = ss.str();
        c->error = mapVectorFieldToError[i];
        c->indices = mapCurvesToClusters.at(i);
        c->curveErrors = mapCurveErrorsToClusters.at(i);
        currentCluster->children.push_back(c);

        c->vectorField = vectorFields.at(i);

        float maxE = -1000;
        for(int j = 0 ; j < (int)c->curveErrors.size() ; ++j){
            float currentError = c->curveErrors.at(j);
            if(currentError > maxE)
                maxE = currentError;
        }
        c->maxError = maxE;
    }

    //
    saveExperiment(outputDirectory, filename, rootCluster);

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
