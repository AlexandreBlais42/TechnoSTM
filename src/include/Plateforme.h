#ifndef PLATEFORME_H
#define PLATEFORME_H

#include <iostream>
#include <string>
#include <cstdint>
#include <cassert>

#include "Vector3D.h"
#include "Serial.h"
#include "Utils.h"

/** @brief Class qui représente la plateforme et permet de la bouger physiquement
 */
class Plateforme{
public:
  Vector3D<uint16_t> position;
  Serial ser;

  Plateforme(const std::string devicePath);

  /** @brief Mets la position de la plateforme relativement à sa position actuelle
   *  @param x La quantité à bouger en x
   *  @param y La quantité à bouger en y
   *  @param z La quantité à bouger en z
   */
  inline void goToRelative(const uint16_t x, const uint16_t y, const uint16_t z);
  /** @brief Mets la position de la plateforme relativement à sa position actuelle
   *  @param coordinates Vector3D spécifiant la quantité à bouger
   */
  inline void goToRelative(const Vector3D<uint16_t> &coordinates);

  /** @brief Mets la position de la plateforme à la position spécifiée
   *  @param x La position en x
   *  @param y La position en y
   *  @param z La position en z
   */
  inline void goToAbsolute(const uint16_t x, const uint16_t y, const uint16_t z);
  /** @brief Mets la position de la plateforme à la position spécifiée
   *  @param coordinates Vector3D spécifiant la position à aller
   */
  inline void goToAbsolute(const Vector3D<uint16_t> &coordinates);

  /** @brief Bouges physiquement la plateforme à la coordonnée dans la variable membre position en l'envoyant sur le port sériel
   */
  void goToPosition() const ;
};

#endif //PLATEFORME_H

