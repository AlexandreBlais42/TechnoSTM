#ifndef GPIO_H
#define GPIO_H

#include <cassert>
#include <cstdint>
#include <iostream>

#ifdef __arm__
#include <wiringPi.h>
#else
typedef enum pinMode_t {
  INPUT,
  OUTPUT,
  PWM_OUTPUT,
  GPIO_CLOCK,
} pinMode_t;
#endif

/** @brief Classe qui gère le GPIO du raspberry pi
 */
class GPIO {
public:
  GPIO();

  /** @brief Mets une pin au mode spécifié
   *  @param pinNumber Le numéro de la pin GPIO
   *  @param mode Le mode de la pin (modes spécifiés dans l'enum pinMode_t)
   */
  static void pinMode(const uint8_t pinNumber, const uint8_t mode);

  /** @brief Set une pin du GPIO du raspberry pi à state
   *  @param pinNumber Le numéro de pin de GPIO
   *  @param state État booléen de la pin à mettre
   */
  static void setPin(const uint8_t pinNumber, const bool state);
};

#endif // GPIO_H
