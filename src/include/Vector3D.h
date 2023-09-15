#ifndef VECTOR3D_H
#define VECTOR3D_H 1

template <typename T>
class Vector3D{
  T x, y, z; 
  Vector3D() : x(0), y(0), z(0){}
  Vector3D(const T x, const T y, const T z) : x(x), y(y), z(z){}
};

#endif //VECTOR3D_H
