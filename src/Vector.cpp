#include "Vector.h"

Vector Vector::operator+(const Vector& rhs) {
    Vector result = *this;
    result += rhs;
    return result;
}

Vector Vector::operator-(const Vector& rhs) {
    Vector result = *this;
    result -= rhs;
    return result;
}

Vector Vector::operator-() {
    Vector result;
    result.x_ = -x_;
    result.y_ = -y_;
    result.z_ = -z_;
    return result;
}

float Vector::normalize() {
    //Returns the distance. Sometimes can be handy.
    float factor = sqrt(
                    x_ * x_
                  + y_ * y_
                  + z_ * z_
                  );

    x_ /= factor;
    y_ /= factor;
    z_ /= factor;

    return factor;
}

Vector Vector::cross(const Vector& rhs) {
    Vector result;
    result.set(y_ * rhs.z_ - z_ * rhs.y_,
               z_ * rhs.x_ - x_ * rhs.z_,
               x_ * rhs.y_ - y_ * rhs.x_);

    return result;
}
