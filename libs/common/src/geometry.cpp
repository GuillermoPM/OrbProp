#include "../include/geometry.h"
#include <math.h>

extern "C"
{
#include "../../cspice/include/SpiceUsr.h"
}
namespace orb{

    Vector3D::Vector3D(){};
    Vector3D::Vector3D(double i, double j, double k) : i(i), j(j), k(k){}
    Vector3D::~Vector3D() = default;

    double Vector3D::norm() const
    {
        return std::sqrt(std::pow(i,2) + std::pow(j,2) + std::pow(k,2));
    }

    Vector3D Vector3D::unit_vector() const
    {
        double magnitude = norm();
        // if (magnitude == 0) {
        //     throw std::runtime_error("Zero vector cannot be normalized");
        // }
        Vector3D uv;
        uv.i = i / magnitude;
        uv.j = j / magnitude;
        uv.k = k / magnitude;

        return uv;
    }

    Vector3D Vector3D::cross(const Vector3D& v1, const Vector3D& v2) {
        Vector3D crossproduct;
        crossproduct.i = v1.j * v2.k - v1.k * v2.j;
        crossproduct.j = v1.k * v2.i - v1.i * v2.k;
        crossproduct.k = v1.i * v2.j - v1.j * v2.i;
        return crossproduct; 
    }

    double Vector3D::dot(const Vector3D& v1, const Vector3D& v2) {
        return v1.i * v2.i + v1.j * v2.j + v1.k * v2.k;
    }

    Vector3D Vector3D::operator+(const Vector3D &other) const
    {
        return Vector3D(i + other.i, j + other.j, k + other.k);
    }

    Vector3D Vector3D::operator-(const Vector3D &other) const
    {
        return Vector3D(i - other.i, j - other.j, k - other.k);
    }

    Vector3D Vector3D::operator*(double scalar) const
    {
        return Vector3D(i * scalar, j * scalar, k * scalar);
    
    }

    Vector3D Vector3D::operator/(double scalar) const
    {
        return Vector3D(i / scalar, j / scalar, k / scalar);
    }

    Vector3D operator*(double scalar, const Vector3D &v)
    {
        return v * scalar;
    }

    Vector3D &Vector3D::operator*=(double scalar)
    {
        i *= scalar;
        j *= scalar;
        k *= scalar;
        return *this;
    }

    Vector3D &Vector3D::operator/=(double scalar)
    {
        i /= scalar;
        j /= scalar;
        k /= scalar;
        return *this;
    }
    Vector3D &Vector3D::operator+=(const Vector3D &other)
    {
        i += other.i;
        j += other.j;
        k += other.k;
        return *this;
    }


}