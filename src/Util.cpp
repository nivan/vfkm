#include "Util.h"
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <queue>
#include "Vector2D.h"

#define xDEBUG

using namespace std;

////////////////////////////////
void Cluster::clearChildren(){
    queue<Cluster*> toProcess;

    int numberOfChildren = children.size();

    for(int i = 0 ; i < numberOfChildren ; ++i){
        toProcess.push(children.at(i));
    }

    while(!toProcess.empty()){
        Cluster* cluster = toProcess.front();
        toProcess.pop();

        delete cluster->vectorField.first;
        delete cluster->vectorField.second;

        int numberOfChildren = cluster->children.size();

        for(int i = 0 ; i < numberOfChildren ; ++i){
            toProcess.push(cluster->children.at(i));
        }

        delete cluster;
    }

    children.clear();
}


////////////////////////////////

Util::Util()
{
}

void Util::loadCurves(string filename, vector<PolygonalPath>& curves,
                      float &xmin, float &xmax, float &ymin, float &ymax,
                      float &tmin, float &tmax){
    cout << "CALLING THIS" << endl;
    xmin = FLT_MAX;
    ymin = FLT_MAX;
    tmin = FLT_MAX;
    xmax = -FLT_MAX;
    ymax = -FLT_MAX;
    tmax = -FLT_MAX;

    // int stride = 0;
    ifstream file (filename.c_str());
    ofstream real_indices("/tmp/real_indices.txt");
    int real_index = 0;
    if (file.is_open())
    {
        //file >> xmin >> xmax >> ymin >> ymax >> tmin >> tmax;
        file >> ymin >> ymax >> xmin >> xmax >> tmin >> tmax;

        vector<pair<Vector2D,float> > curveContents;

        while (file.good()  && !file.eof())
        {
            float x,y,t;
            //file >> x >> y >> t;
            file >> y >> x >> t;

            if(x == 0 && y == 0 && t == 0) {
                if(curveContents.size() >=2){
                    real_indices << real_index << endl;
                    curves.push_back(PolygonalPath(curveContents));
                }
                real_index++;
                curveContents.clear();
            } else if (x < xmin || x > xmax ||
                  y < ymin || y > ymax ||
                  t < tmin || t > tmax) {//end of a curve
                if(curveContents.size() >=2){
                    real_indices << real_index << endl;
                    curves.push_back(PolygonalPath(curveContents));
                }
                curveContents.clear();
            } else {
                pair<Vector2D, float> newPoint = make_pair(Vector2D(x,y), t);
                if (curveContents.size() == 0)
                    curveContents.push_back(newPoint);
                else if (t == curveContents.back().second)
                    continue;
                else if (x == curveContents.back().first.X() &&
                    y == curveContents.back().first.Y())
                    continue;
                else
                    curveContents.push_back(newPoint);
            }
        }
        file.close();
    }
    else{
        cerr << "Unable to open file " << filename << endl;
    }


    int numberOfCurvesRead = curves.size();
    ofstream outfile("msr_campus_cropped.txt");
    outfile << xmin << " " << xmax << " " << ymin << " " << ymax << " " << tmin << " " << tmax;
    for(int i = 0 ; i < numberOfCurvesRead ; ++i){
        //cout << "Curve " << i << " = " << curves.at(i).toString() << endl;
        PolygonalPath &curveContents = curves.at(i);
        int curveSize = curveContents.numberOfPoints();

        for(int j = 0 ; j < curveSize ; ++j){
            pair<Vector2D,float> pointTime = curveContents.getPoint(j);
            outfile << pointTime.first.X() << " " << pointTime.first.Y() << " " << pointTime.second << endl;
        }
        outfile << "0 0 0" << endl;
    }





#ifdef DEBUG
    int numberOfCurvesRead = curves.size();

    cout << "numberOfCurvesRead = " << numberOfCurvesRead << endl;

    for(int i = 0 ; i < numberOfCurvesRead ; ++i){
        cout << "Curve " << i << " = " << curves.at(i).toString() << endl;
    }
#endif   
}

void Util::loadCurves(std::string filename, std::vector<PolygonalPath>& curves){
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float tmin;
    float tmax;

    loadCurves(filename, curves, xmin,xmax,ymin,ymax,tmin,tmax);
}

void lat_long_to_mercator(float &x, float &y,
                          float latitude,
                          float longitude)
{
    float lat_rad = latitude * (M_PI / 180.f);
    y = logf(1.0/cosf(lat_rad) + tanf(lat_rad));
    x = longitude * (M_PI / 180.f);
}

#include <cstdlib>

void Util::loadCurvesAndProject(std::string filename, std::vector<PolygonalPath>& curves,
                                 float &xmin, float &xmax, float &ymin, float &ymax,
                                 float &tmin, float &tmax){
    xmin = FLT_MAX;
    ymin = FLT_MAX;
    tmin = FLT_MAX;
    xmax = -FLT_MAX;
    ymax = -FLT_MAX;
    tmax = -FLT_MAX;

     cout << "LOADING " << filename << endl;

    ifstream file (filename.c_str());
    if (file.is_open())
    {
        string line;
        vector<pair<Vector2D,float> > curveContents;

        while ( file.good() && !file.eof())
        {
            float x,y,t;
            int count=0;
            file >> x >> y >>t;

            if(x == 0 && y == 0 && t ==0){//end of a curve
                if(curveContents.size() >=2){
                    curves.push_back(PolygonalPath(curveContents));
                }

                curveContents.clear();
            }
            else{
                if (count++ % 10 != 0)
                    continue;
                float xx = x;
                float yy = y;
                // float tt = t;

                lat_long_to_mercator(x,y,xx,yy);

                pair<Vector2D, float> newPoint = make_pair(Vector2D(x,y), t);
                curveContents.push_back(newPoint);


                if(xmin > x)
                    xmin = x;
                if(xmax < x)
                    xmax = x;
                if(ymin > y)
                    ymin = y;
                if(ymax < y)
                    ymax = y;
                if(tmin > t)
                    tmin = t;
                if(tmax < t)
                    tmax = t;
            }
        }
        file.close();
    }
    else{
        cerr << "Unable to open file " << filename << endl;
        exit(1);
    }


#ifdef DEBUG
    int numberOfCurvesRead = curves.size();

    cout << "numberOfCurvesRead = " << numberOfCurvesRead << endl;

    for(int i = 0 ; i < numberOfCurvesRead ; ++i){
        cout << "Curve " << i << " = " << curves.at(i).toString() << endl;
    }
#endif
}

#define rads_over_degrees (3.1415926535897931 / 180.0)
#include <cmath>

void Util::to_mercator(const float &lat, const float &lon, float &xMerc, float &yMerc){
    float latitude = lat * rads_over_degrees;
    xMerc = lon * rads_over_degrees;
    yMerc = log(tan(latitude) + 1.0/cos(latitude));
}
