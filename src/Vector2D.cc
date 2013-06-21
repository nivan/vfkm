#include "Vector2D.h"
#include <sstream>
#include <cstdlib>
#include <iostream>

using namespace std;

#ifdef DEBUG
Vector2D::Vector2D(float x, float y):
    x(x),
    y(y)
{}

float Vector2D::X() const{
    return x;
}

float Vector2D::Y() const{
    return y;
}

double Vector2D::dot(float vx, float vy) const{
    return x * vx + y * vy;
}

double Vector2D::dot(const Vector2D &v) const{
    return x * v.X() + y * v.Y();
}

float Vector2D::length() const{
    return sqrt(dot(*this));
}

double Vector2D::length2() const{
    return dot(*this);
}

void Vector2D::normalize(){
    float l = length();

    if(l > 0)
	scale(1.0f/l);
}

void Vector2D::scale(float lambda){
    x *= lambda;
    y *= lambda;
}

void Vector2D::add(float vx, float vy){
    x += vx;
    y += vy;
}

void Vector2D::add(const Vector2D &v){
    x += v.X();
    y += v.Y();
}

void Vector2D::subtract(float vx, float vy){
    x -= vx;
    y -= vy;
}

void Vector2D::subtract(const Vector2D &v){
    x -= v.X();
    y -= v.Y();
}

bool Vector2D::equals(const Vector2D &v){
    return ((x == v.X())&&(y == v.Y()));
}

Vector2D Vector2D::operator-(Vector2D v){
    return Vector2D(this->x - v.X(), this->y - v.Y());
}

Vector2D Vector2D::operator+(Vector2D v){
    return Vector2D(this->x + v.X(), this->y + v.Y());
}

Vector2D Vector2D::operator* (float lambda){
    return Vector2D(lambda * x, lambda * y);
}
#endif


void Vector2D::rotate(float angle, RotationOrientation orientation){
    
    float coef = (orientation == Vector2D::COUNTER_CLOCK_WISE)?-1.0:1.0;

    float c = cos(angle);
    float s = sin(angle);

    float newX = c * x + coef * s * y;
    float newY = (-coef) * s * x + c * y;

    x = newX;
    y = newY;
}

bool Vector2D::compareLowerLeft(const Vector2D& v1, const Vector2D& v2){
    return (v1.X() < v2.X()) || ((v1.X() == v2.X()) && (v1.Y() < v2.Y()));
}

Vector2D::RelativePosition Vector2D::isToTheLeft(const Vector2D& v1, const Vector2D& v2){
    float zCoord =  (v2.X() * v1.Y() - v2.Y() * v1.X());
    
    if (zCoord == 0.0)
	return Vector2D::ALIGNED;
    else if(zCoord > 0.0)
	return Vector2D::TO_THE_LEFT;
    else
	return Vector2D::TO_THE_RIGHT;
}

Vector2D::RelativePosition Vector2D::isALeftTurn(const Vector2D& v1, const Vector2D& v2, const Vector2D& v3){
    Vector2D v12(v2.X()-v1.X(), v2.Y()-v1.Y());
    Vector2D v23(v3.X()-v2.X(), v3.Y()-v2.Y());
    return isToTheLeft(v23,v12);
}

string Vector2D::toString(int precision) const{
    stringstream ss;
    ss.precision(precision);
    ss << "(" << x << "," << y << ")";

    return ss.str();
}

string Vector2D::toStringRelativePosition(RelativePosition rp){
    switch(rp){
    case(Vector2D::TO_THE_LEFT):
	return "TO_THE_LEFT";
	break;
    case(Vector2D::TO_THE_RIGHT):
	return "TO_THE_RIGHT";
	break;
    case(Vector2D::ALIGNED):
	return "ALIGNED";
	break;
    default:
	cerr << "Error in string Vector2D::toStringRelativePosition(RelativePosition rp)" << endl;
	exit(1);
    };
}


std::string Vector2D::toString(RelativePosition rp){
    switch(rp){
    case(Vector2D::TO_THE_LEFT):
        return string("to the left");
        break;
    case(Vector2D::TO_THE_RIGHT):
        return string("to the right");
        break;
    case(Vector2D::ALIGNED):
        return string("aligned");
        break;
    default:
        cerr << "Error invalid realtive position" << endl;
        exit(1);
    };
}

float Vector2D::dot(Vector2D v, Vector2D u){
    return ((v.X() * u.X()) + (v.Y() * u.Y()));
}

bool Vector2D::intersect(const Vector2D &p1, const Vector2D &p2,
                         const Vector2D &q1, const Vector2D &q2,
                         float &lambda1, float &lambda2){

    Vector2D::RelativePosition rpP1 = Vector2D::isALeftTurn(p1, p2, q1);
    Vector2D::RelativePosition rpP2 = Vector2D::isALeftTurn(p1,p2,q2);

    bool differentSides1 = ((rpP1 == Vector2D::TO_THE_LEFT) && (rpP2 == Vector2D::TO_THE_RIGHT)) ||
                           ((rpP1 == Vector2D::TO_THE_RIGHT) && (rpP2 == Vector2D::TO_THE_LEFT));

    Vector2D::RelativePosition rpQ1 = Vector2D::isALeftTurn(q1, q2, p1);
    Vector2D::RelativePosition rpQ2 = Vector2D::isALeftTurn(q1, q2, p2);

    bool differentSides2 = ((rpQ1 == Vector2D::TO_THE_LEFT) && (rpQ2 == Vector2D::TO_THE_RIGHT)) ||
                           ((rpQ1 == Vector2D::TO_THE_RIGHT) && (rpQ2 == Vector2D::TO_THE_LEFT));


    if(differentSides1 && differentSides2){
        //the segments intersect
        //find lambdas

        double det =  (p2.X() - p1.X()) * (q1.Y() - q2.Y()) - (p2.Y() - p1.Y()) * (q1.X() - q2.X());
        double detx = (q1.X() - p1.X()) * (q1.Y() - q2.Y()) - (q1.Y() - p1.Y()) * (q1.X() - q2.X());
        double dety = (p2.X() - p1.X()) * (q1.Y() - p1.Y()) - (p2.Y() - p1.Y()) * (q1.X() - p1.X());

        lambda1 = detx / det;
        lambda2 = dety / det;

        return true;
    }
    else{
        lambda1 = -1;
        lambda2 = -1;
        return false;
    }

}
