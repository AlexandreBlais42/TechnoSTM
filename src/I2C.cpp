#include "I2C.h"

I2C::I2C(const uint8_t deviceAddr) : deviceAddr(deviceAddr) {
#ifdef __arm__
  file = std::fstream("/dev/i2c-" + std::to_string(deviceAddr), file.binary | file.trunc | file.in | file.out);
  if (!file.is_open()){
    std::cout << "Il n'y a pas d'appareil à l'addresse " << std::to_string(deviceAddr) << std::endl;
    exit(1);
  }
#else
  std::cout << "I2C::I2C appellé avec deviceAddr = " << std::to_string(deviceAddr) << std::endl;
#endif // __arm__
}

I2C::~I2C(){
#ifdef __arm__
  file.close();
#else
  std::cout << "I2C::~I2C appellé " << std::endl;
#endif // __arm__
}

