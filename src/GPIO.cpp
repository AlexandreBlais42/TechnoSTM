#include "GPIO.h"

GPIO::GPIO(){}

void GPIO::setPin(const uint8_t pinNumber, const bool state){
  // À implémenter pour interfacer avec le raspberry pi
  assert(pinNumber >=2 && pinNumber <= 27);
  std::cout << "GPIO::setPin appellé avec pinNumber = " << pinNumber << " state = " << state << std::endl;
}
