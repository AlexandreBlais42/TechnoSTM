#ifndef STEPMOTOR_H
#define STEPMOTOR_H

#include <cstdint>

#include "Utils.h"
#include "GPIO.h"

class StepMotor{
public:
  int32_t position;
  int8_t steps;

  StepMotor(const int32_t position);

  /** @brief Mets la position de la plateforme relativement à sa position actuelle 
   *  @param position La quantité à bouger
   */
  void goToRelative(const int32_t pos);
  /** @brief Mets la position de la plateforme à la position spécifiée
   *  @param position La position à aller
   */
  inline void goToAbsolute(const int32_t pos);
};

#endif //STEPMOTOR_H
