#include "ADS1115.h"

ADS1115::ADS1115(const uint8_t deviceAddr)
    : I2C(deviceAddr), configs(DEFAULT_ADS1115_CONFIGS),
      addressPointer(Conversion) {}


void ADS1115::setOperationalStatus(const bool status){
  configs &= OPERATIONAL_STATUS_MASK;
  configs |= ((uint16_t) status) << OPERATIONAL_STATUS_OFFSET;
}

void ADS1115::setChannel(const uint8_t channel) {
  assert(channel < 4);
  configs &= CHANNEL_MASK;
  configs |= ((uint16_t) channel + 0b100) << CHANNEL_OFFSET;
}

void ADS1115::setProgrammableGain(const uint8_t gain){
  assert(gain < 8);
  configs &= PROGRAMMABLE_GAIN_MASK;
  configs |= ((uint16_t) gain) << PROGRAMMABLE_GAIN_OFFSET;
}

void ADS1115::setSingleConversionMode(const bool mode){
  configs &= SINGLE_CONVERSION_MASK;
  configs |= ((uint16_t) mode) << SINGLE_CONVERSION_OFFSET;
}

void ADS1115::setDataRate(const uint8_t rate){
  assert(rate < 8);
  configs &= DATA_RATE_MASK;
  configs |= ((uint16_t) rate) << DATA_RATE_OFFSET;
}

void ADS1115::setCompareMode(const bool mode){
  configs &= COMPARE_MODE_MASK;
  configs |= ((uint16_t) mode) << COMPARE_MODE_OFFSET;
}

void ADS1115::setComparatorPolarity(const bool mode){
  configs &= COMPARATOR_POLARITY_MASK;
  configs |= ((uint16_t) mode) << COMPARATOR_POLARITY_OFFSET;
}

void ADS1115::setLatchingComparator(const bool mode){
  configs &= LATCHING_COMPARATOR_MASK;
  configs |= ((uint16_t) mode) << LATCHING_COMPARATOR_OFFSET;
}

void ADS1115::setComparatorQueue(const uint8_t mode){
  assert(mode < 4);
  configs &= COMPARATOR_QUEUE_MASK;
  configs |= ((uint16_t) mode) << COMPARATOR_QUEUE_OFFSET;
}

void ADS1115::writeConfigs(){
  writeReg16(Config, configs); 
  addressPointer = Config;
}

uint16_t ADS1115::read() {
  if (addressPointer != Conversion){
    writeData<uint8_t>(Conversion);
    addressPointer = Conversion;
  }
  return readData<uint16_t>();
}
