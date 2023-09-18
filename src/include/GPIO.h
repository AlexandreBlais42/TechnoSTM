#ifndef GPIO_H
#define GPIO_H

#include <iostream>
#include <cstdint>
#include <cassert>

class GPIO{
public:
  GPIO();

  /** @brief Set une pin du GPIO du raspberry pi à state
   *  @param pinNumber Le numéro de pin de GPIO
   *  @param state État booléen de la pin à mettre
   */
  static void setPin(const uint8_t pinNumber, const bool state);
};

#endif //GPIO_H
