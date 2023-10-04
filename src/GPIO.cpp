#include "GPIO.h"
#include <string>

GPIO::GPIO() {}

void GPIO::pinMode(const uint8_t pinNumber, const uint8_t mode) {
#ifdef __arm__
  pinMode(pinNumber, mode);
#else
  std::cout << "GPIO::setPin appellé avec pinNumber = "
            << std::to_string(pinNumber) << " mode = " << mode << std::endl;
#endif // __arm__
}

void GPIO::setPin(const uint8_t pinNumber, const bool state) {
  assert(pinNumber >= 2 && pinNumber <= 27);
#ifdef __arm__
  digitalWrite(pinNumber, state);
#else
  std::cout << "GPIO::setPin appellé avec pinNumber = "
            << std::to_string(pinNumber) << " state = " << state << std::endl;
#endif // __arm__
}
