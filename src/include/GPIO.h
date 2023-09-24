#ifndef GPIO_H
#define GPIO_H

#include <cassert>
#include <cstdint>
#include <iostream>

#ifdef __arm__
#include <wiringPi.h>
#endif // __arm__

/** @brief Classe qui gère le GPIO du raspberry pi
 */
class GPIO {
public:
  GPIO();

  /** @brief Set une pin du GPIO du raspberry pi à state
   *  @param pinNumber Le numéro de pin de GPIO
   *  @param state État booléen de la pin à mettre
   */
  static void setPin(const uint8_t pinNumber, const bool state);
};

#endif // GPIO_H
