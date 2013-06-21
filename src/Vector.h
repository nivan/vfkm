#pragma once

#include <string>
#include <cmath>

#define VECTOR_TYPE float

class Vector
{
private:
    int dimension;
    VECTOR_TYPE *values;
public:
    Vector() { dimension=0; values = NULL; }
    explicit Vector(int dimension);
    Vector(const Vector& v);
    ~Vector();

    inline int getDimension() const { return dimension; }
    inline void setValue(int index, VECTOR_TYPE value) {
        values[index] = value;
    }
    inline void setValues(VECTOR_TYPE *newValues) {
      for(int i = 0 ; i < dimension ; ++i)
        values[i] = newValues[i];
    }
    void alloc(int size) {
        if (values) {
            delete[] values;
        }
        values = new VECTOR_TYPE[size];
        dimension = size;
    }
    inline void setValues(const Vector& v) {
      for(int i = 0 ; i < dimension ; ++i)
        values[i] = v[i];
    }
    inline void setValues(VECTOR_TYPE newValue) {
      for(int i = 0 ; i < dimension ; ++i)
        values[i] = newValue;
    }

    inline VECTOR_TYPE operator[](int index) const {
        return values[index];
    }

    inline VECTOR_TYPE& operator[](int index) {
        return values[index];
    }

    inline Vector& operator=(const Vector& v) { 
        if (v.getDimension() != getDimension()) {
            alloc(v.getDimension());
        }
        setValues(v);
        return *this;
    }

    inline void add(const Vector& v) {
      for(int i = 0 ; i < dimension ; ++i)
        values[i] += v[i];
    }
    inline void add(const Vector& v, Vector& result) {
      for(int i = 0 ; i < dimension ; ++i)
        result.setValue(i,values[i] + v[i]);
    }
    inline void scale(VECTOR_TYPE s) {
      for(int i = 0 ; i < dimension ; ++i)
	values[i] *= s;
    }
    inline void scale(VECTOR_TYPE s, Vector& v) {
      for(int i = 0 ; i < dimension ; ++i)
        v.setValue(i, values[i] * s);
    }
    inline VECTOR_TYPE length2() const {
      VECTOR_TYPE norm = 0.0;
      for(int i = 0 ; i < dimension ; ++i)
        norm += (values[i] * values[i]);
      return norm;
    }
    inline VECTOR_TYPE length() const { return sqrtf(length2()); }

    inline VECTOR_TYPE dot(const Vector& v) const {
      VECTOR_TYPE norm = 0.0;
      for(int i = 0 ; i < dimension ; ++i)
        norm += (values[i] * v[i]);
      return norm;
    }

    inline void add_scale(const Vector &v, VECTOR_TYPE b, VECTOR_TYPE a=1.0) {
        for (int i=0; i<dimension; ++i) {
            values[i] = a * values[i] + b * v[i];
        }
    }

    inline Vector &operator/=(const Vector &other) {
        for (int i=0; i<dimension; ++i) {
            values[i] /= other.values[i];
        }
        return *this;
    }

    inline void sub(const Vector& v) {
      for(int i = 0 ; i < dimension ; ++i)
        values[i] -= v[i];
    }
    inline void sub(const Vector& v, Vector& result) {
      for(int i = 0 ; i < dimension ; ++i)
        result.setValue(i,values[i] - v[i]);
    }

    inline Vector &operator+=(const Vector &v) {
        add(v);
        return *this;
    }

    inline Vector &operator-=(const Vector &v) {
        sub(v);
        return *this;
    }

    std::string toString() const;
};

inline Vector operator-(const Vector &v1, const Vector &v2)
{
    Vector result(v1);
    return result -= v2;
}

inline Vector operator+(const Vector &v1, const Vector &v2)
{
    Vector result(v1);
    return result += v2;
}

inline Vector operator*(VECTOR_TYPE k, const Vector &v)
{
    Vector result(v);
    result.scale(k);
    return result;
}

inline Vector operator*(const Vector &v, VECTOR_TYPE k)
{
    Vector result(v);
    result.scale(k);
    return result;
}

inline Vector operator/(const Vector &v1, const Vector &v2)
{
    Vector result(v1);
    result /= v2;
    return result;
}
