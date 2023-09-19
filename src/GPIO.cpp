#include "GPIO.h"

GPIO::GPIO(){}

void GPIO::setPin(const uint8_t pinNumber, const bool state){
  assert(pinNumber >=2 && pinNumber <= 27);
#ifdef __arm__
  // À implémenter pour interfacer avec le raspberry pi
#else
  std::cout << "GPIO::setPin appellé avec pinNumber = " << pinNumber << " state = " << state << std::endl;
#endif __arm__
}
