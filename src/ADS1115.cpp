#include "ADS1115.h"

ADS1115::ADS1115(const uint8_t deviceAddr)
    : I2C(deviceAddr), configs(DEFAULT_ADS1115_CONFIGS),
      addressPointer(Conversion) {}

void ADS1115::setChannel(const uint8_t channel) {
  assert(channel <= 3);
  configs &= CHANNEL_MASK;
  configs |= ((uint16_t) channel + 0b100) << 12;
}

uint16_t ADS1115::read() {
#ifdef __arm__
  return I2C::read<uint16_t>();
#else
  std::cout << "ADS1115::read() appellé" << std::endl;
  return 0;
#endif // __arm__
}
