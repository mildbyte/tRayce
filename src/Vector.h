//A 3D vector, used for storing colors and coordinates; vector algebra.
#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>

class Vector {
private:
    double x_, y_, z_;
public:
    //Constructor
    Vector(){}
    Vector(double x, double y, float z) {x_ = x; y_ = y; z_ = z;}

    //Getters
    double getX() {return x_;}
    double getY() {return y_;}
    double getZ() {return z_;}

    //Setters
    void setX(double x) {x_ = x;}
    void setY(double y) {y_ = y;}
    void setZ(double z) {z_ = z;}
    void set(double x, double y, float z) {x_ = x; y_ = y; z_ = z;}

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

    //Dot and cross products
    double dot(const Vector& rhs) {return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;}
    Vector cross(const Vector& rhs); //too long to be featured here

    double normalize();
    double modulus() {return sqrtf(x_ * x_ + y_ * y_ + z_ * z_);}

    //Reinhard mapping (0..inf to 0..1)
    void reinhardMap() {
        x_ = x_ / (1 + x_);
        y_ = y_ / (1 + y_);
        z_ = z_ / (1 + z_);
    }

    void toPow(double pow) {
        x_ = powf(x_, pow);
        y_ = powf(y_, pow);
        z_ = powf(z_, pow);
    }

};

#endif
