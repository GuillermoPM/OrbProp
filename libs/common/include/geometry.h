#pragma once
namespace orb
{
    
    class Vector3D
    {
    public:
        double i, j, k;

        Vector3D();
        Vector3D(double i, double j, double k);
        ~Vector3D();

        double norm() const;
        Vector3D unit_vector() const;
        static Vector3D cross(const Vector3D &v1, const Vector3D &v2);
        static double dot(const Vector3D &v1, const Vector3D &v2);

        Vector3D operator+(const Vector3D &other) const;
        Vector3D operator-(const Vector3D &other) const;
        Vector3D operator*(double scalar) const;
        Vector3D operator/(double scalar) const;
        Vector3D &operator*=(double scalar);
    };
    Vector3D operator*(double scalar, const Vector3D &v);


} // namespace orb
