#ifndef GPIO_H
#define GPIO_H

#include <iostream>
#include <cstdint>

class GPIO{
public:
  GPIO();

  static void setPin(const uint8_t pinNumber, const bool state);
};

#endif //GPIO_H
