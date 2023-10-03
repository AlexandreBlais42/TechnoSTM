#ifndef STEPMOTOR_H
#define STEPMOTOR_H

#include <cstdint>
#include <array>

#include "GPIO.h"
#include "Utils.h"

/** @brief Classe qui gère les mouvements du step moteur
 */
class StepMotor : private GPIO {
public:
  int32_t position;
  std::array<uint8_t, 4> pins;
  int8_t pinIndex;

  /** @brief Constructeur de la classe StepMotor
   *  @param position La position initiale du moteur
   *  @param pins Une liste de 4 valeurs pour les pins
   */
  StepMotor(const int32_t position, const std::array<uint8_t, 4> pins);

  /** @brief Mets la position de la plateforme relativement à sa position
   * actuelle
   *  @param pos La quantité à bouger
   */
  void goToRelative(const int32_t pos);
  /** @brief Mets la position de la plateforme à la position spécifiée
   *  @param pos La position à aller
   */
  inline void goToAbsolute(const int32_t pos);
};

#endif // STEPMOTOR_H
