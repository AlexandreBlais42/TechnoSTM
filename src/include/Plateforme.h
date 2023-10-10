#ifndef PLATEFORME_H
#define PLATEFORME_H

#include <cstdint>
#include <iostream>
#include <string>

#include "Serial.h"
#include "Utils.h"
#include "Vector3D.h"

/** @brief Class qui représente la plateforme et permet de la bouger
 * physiquement
 */
class Plateforme : private Serial {
public:
  Vector3D<uint16_t> position;

  Plateforme(const std::string devicePath);

  /** @brief Bouge plateforme relativement à sa position
   * actuelle
   *  @param x La quantité à bouger en x
   *  @param y La quantité à bouger en y
   *  @param z La quantité à bouger en z
   */
  inline void setPositionRelative(const uint16_t x, const uint16_t y,
                                  const uint16_t z);
  /** @brief Bouge plateforme relativement à sa position
   * actuelle
   *  @param coordinates Vector3D spécifiant la quantité à bouger
   */
  inline void setPositionRelative(const Vector3D<uint16_t> &coordinates);

  /** @brief Bouge la plateforme à la position spécifiée
   *  @param x La position en x
   *  @param y La position en y
   *  @param z La position en z
   */
  inline void setPositionAbsolute(const uint16_t x, const uint16_t y,
                                  const uint16_t z);
  /** @brief Bouge plateforme à la position spécifiée
   *  @param coordinates Vector3D spécifiant la position à aller
   */
  inline void setPositionAbsolute(const Vector3D<uint16_t> &coordinates);

  /** @brief Bouges physiquement la plateforme à la coordonnée dans la variable
   * membre position en envoyant des commandes sur le port sériel
   */
  void moveToPosition() const;
};

#endif // PLATEFORME_H
