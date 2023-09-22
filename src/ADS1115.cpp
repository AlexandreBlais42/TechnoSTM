#include "ADS1115.h"


ADS1115::ADS1115(){

}

uint16_t ADS1115::read(const uint8_t channel) const {
  assert(channel <= 3);

#ifdef __arm__
  // À implémenter
  return 0;
#else
  std::cout << "ADS1115::read() appellé avec channel = " << (uint16_t) channel << std::endl;
  return 0;
#endif // __arm__
}
