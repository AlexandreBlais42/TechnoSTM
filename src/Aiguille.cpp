#include "Aiguille.h"

Aiguille::Aiguille(const uint8_t deviceAddr) : ADS1115(deviceAddr) {}

#ifdef __arm__
int16_t Aiguille::readVoltage() { return read(); }
#else
int16_t Aiguille::readVoltage(){
  int16_t toReturn;
  std::cin >> toReturn;
  return toReturn;
}
#endif // __arm__
