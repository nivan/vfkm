#ifndef VECTOR_2D
#define VECTOR_2D

#include <cmath>
#include <string>

class Vector2D {
public:
    float x;
    float y;

public:
    enum RotationOrientation{
	CLOCK_WISE,
	COUNTER_CLOCK_WISE
    };

    enum RelativePosition {
	TO_THE_LEFT,
	TO_THE_RIGHT,
	ALIGNED
    };

public:
#ifdef DEBUG
    Vector2D(float x = 0.0f, float y = 0.0f);
    float X() const;
    float Y() const;
    double dot(float vx, float vy) const;
    double dot(const Vector2D &v) const;
    float length() const;
    double length2() const;
    void normalize();
    void scale(float lambda);
    void add(float x, float y);
    void add(const Vector2D &v);
    void subtract(float x, float y);
    void subtract(const Vector2D &v);
    bool equals(const Vector2D &v);
    Vector2D operator- (Vector2D);
    Vector2D operator+ (Vector2D);
    Vector2D operator* (float lambda);
#else
    inline Vector2D(float x=0.0f, float y=0.0f):x(x), y(y) {}
    inline float X() const { return x; }
    inline float Y() const { return y; }
    inline double dot(float vx, float vy) const {
      return x * vx + y * vy;
    }
    inline double dot(const Vector2D &v) const {
      return x * v.X() + y * v.Y();
    }
    inline float length() const {
      return sqrt(x*x + y*y);
    }
    inline double length2() const {
      return x*x + y*y;
    }
    inline void normalize() {
      float l = length();
      if (l > 0)
	scale(1.0f/l);
    }
    inline void scale(float lambda) {
      x *= lambda;
      y *= lambda;
    }
    inline void add(float vx, float vy) {
      x += vx;
      y += vy;
    }
    inline void add(const Vector2D &v) {
      x += v.X();
      y += v.Y();
    }
    inline void subtract(float vx, float vy) {
      x -= vx;
      y -= vy;
    }
    inline void subtract(const Vector2D &v) {
      x -= v.X();
      y -= v.Y();
    }
    inline bool equals(const Vector2D &v) {
      return ((x==v.X()) && (y==v.Y()));
    }
    inline Vector2D operator-(Vector2D v){
      return Vector2D(this->x - v.X(), this->y - v.Y());
    }
    inline Vector2D operator+(Vector2D v){
      return Vector2D(this->x + v.X(), this->y + v.Y());
    }
    inline Vector2D operator* (float lambda){
      return Vector2D(lambda * x, lambda * y);
    }
    inline float distance(const Vector2D &v) {
        Vector2D x(*this);
        x.subtract(v);
        return x.length();
    }
#endif

    void rotate(float angle, RotationOrientation orientation = COUNTER_CLOCK_WISE);

    std::string toString(int precision = 4) const;
    static std::string toString(RelativePosition);
    static std::string toStringRelativePosition(RelativePosition rp);

    static bool compareLowerLeft(const Vector2D& v1, const Vector2D& v2);    
    static RelativePosition isToTheLeft(const Vector2D& v1, const Vector2D& v2);
    static RelativePosition isALeftTurn(const Vector2D& v1, const Vector2D& v2, const Vector2D& v3);
    static bool intersect(const Vector2D &p1, const Vector2D &p2,
                          const Vector2D &q1, const Vector2D &q2,
                          float &lambda1, float &lambda2);
    static float dot(Vector2D v, Vector2D u);
    static float cross(const Vector2D &v1, const Vector2D &v2) {
        return v1.x * v2.y - v2.x * v1.y;
    }

    static Vector2D lerp(const Vector2D &v1, const Vector2D &v2, float u) {
        Vector2D result(v2);
        result.scale(u);
        Vector2D tmp(v1);
        tmp.scale(1-u);
        result.add(tmp);
        return result;
    }
};

#endif
