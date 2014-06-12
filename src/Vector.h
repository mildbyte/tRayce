//A 3D vector, used for storing colors and coordinates; vector algebra.
#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <cstdio>
#include <stdexcept>

const double PI = 4 * atan(1);
const double EPSILON = 0.0000001;

class Vector {
private:
    double x_, y_, z_;
public:
    //Constructor
    Vector(){}
    Vector(double x, double y, double z) {x_ = x; y_ = y; z_ = z;}

    //Getters
    inline double getX() {return x_;}
    inline double getY() {return y_;}
    inline double getZ() {return z_;}

    //Setters
    inline void setX(double x) {x_ = x;}
    inline void setY(double y) {y_ = y;}
    inline void setZ(double z) {z_ = z;}
    inline void set(double x, double y, double z) {x_ = x; y_ = y; z_ = z;}

    //Algebra
    Vector& operator=(const Vector& rhs) {x_ = rhs.x_; y_ = rhs.y_; z_ = rhs.z_; return *this;}

    //Compound operators
    Vector& operator+=(const Vector& rhs) {x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_; return *this;}
    Vector& operator-=(const Vector& rhs) {x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_; return *this;}
    Vector& operator/=(double rhs) {x_ /= rhs; y_ /= rhs; z_ /= rhs; return *this;}
    Vector& operator*=(double rhs) {x_ *= rhs; y_ *= rhs; z_ *= rhs; return *this;}

    //Binary and unary operators
    Vector operator+(const Vector& rhs);
    Vector operator-(const Vector& rhs);
    Vector operator-();

    Vector operator*(const double rhs) {return Vector(x_ * rhs, y_ * rhs, z_ * rhs);}
    
    double& operator[](int ix) {
        switch (ix) {
            case 0: return x_; break;
            case 1: return y_; break;
            case 2: return z_; break;
            default: throw std::out_of_range("vector subscript");
        }
    }

    //Dot and cross products
    double dot(const Vector& rhs) {return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;}
    Vector cross(const Vector& rhs); //too long to be featured here

    double normalize();
    double modulus() {return sqrt(x_ * x_ + y_ * y_ + z_ * z_);}

    void print() {printf("%f, %f, %f", x_, y_, z_);}

    //Reinhard mapping (0..inf to 0..1)
    void reinhardMap() {
        x_ = x_ / (1 + x_);
        y_ = y_ / (1 + y_);
        z_ = z_ / (1 + z_);
    }

    void toPow(double p) {
        x_ = pow(x_, p);
        y_ = pow(y_, p);
        z_ = pow(z_, p);
    }

};

#endif
