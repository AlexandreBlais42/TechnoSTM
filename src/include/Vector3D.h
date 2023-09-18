#ifndef VECTOR3D_H
#define VECTOR3D_H

template <typename T>
/** @brief Une structure template qui contient x, y, z et qui permet de les additionner facilement
 *
 */
class Vector3D{
public:
  T x, y, z; 
  Vector3D() : x(0), y(0), z(0){}
  Vector3D(const T x, const T y, const T z) : x(x), y(y), z(z){}

  /** @brief Additione un Vector3D à l'objet
   *  @param other L'autre objet 
   */
  void operator+=(const Vector3D<T> &other){
    x += other.x;
    y += other.y;
    z += other.z;
  }

  /** @brief Soustrait un Vector3D à l'objet
   *  @param other L'autre objet
   */
  void operator-=(const Vector3D<T> &other){
    x -= other.x;
    y -= other.y;
    z -= other.z;
  }
};


template <typename T>
/** @brief Additionne deux objets Vector3D
 *  @param a Le premier objet
 *  @param b Le deuxième objet
 *  @return La somme des deux objets
 */
Vector3D<T> operator+(const Vector3D<T> &a, const Vector3D<T> &b){
  return Vector3D<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename T>
/** @brief Soustrait deux objets Vector3D
 *  @param a Le premier objet
 *  @param b Le deuxième objet
 *  @return La différence des deux objets
 */
Vector3D<T> operator-(const Vector3D<T> &a, const Vector3D<T> &b){
  return Vector3D<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}

#endif //VECTOR3D_H
