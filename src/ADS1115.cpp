#include "ADS1115.h"
#define __arm__

ADS1115::ADS1115(const uint8_t deviceAddr) : I2C(deviceAddr), configs(DEFAULT_ADS1115_CONFIGS) {

}

void ADS1115::setChannel(const uint8_t channel){
#ifdef __arm__
  // I2C::write();
#else
  std::cout << "ADS1115::setChannel() appellé avec channel = " << (uint16_t) channel << std::endl;
#endif // __arm__

}

uint16_t ADS1115::read(){
#ifdef __arm__
  return I2C::read<uint16_t>();
#else
  std::cout << "ADS1115::read() appellé" << std::endl;
  return 0;
#endif // __arm__
}
