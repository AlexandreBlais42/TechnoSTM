#include "GPIO.h"

GPIO::GPIO() {}

void GPIO::begin(){
  #ifdef __arm__
  wiringPiSetupGpio();
  #endif // __arm__
}

void GPIO::_pinMode(const uint8_t pinNumber, const uint8_t mode) {
  assert(pinNumber <= 27);
#ifdef __arm__
  pinMode(pinNumber, mode);
#else
  DEBUG << "GPIO::setPin appellé avec pinNumber = "
            << std::to_string(pinNumber) << " mode = " << (mode ? "output" : "Input") << std::endl;
#endif // __arm__
}

void GPIO::setPin(const uint8_t pinNumber, const bool state) {
  assert(pinNumber <= 27);
#ifdef __arm__
  digitalWrite(pinNumber, state);
#else
  DEBUG << "GPIO::setPin appellé avec pinNumber = "
            << std::to_string(pinNumber) << " state = " << state << std::endl;
#endif // __arm__
}
