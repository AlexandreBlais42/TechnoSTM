#ifndef ADS1115_H
#define ADS1115_H

#include <iostream>
#include <cstdint>
#include <cassert>

#include "I2C.h"

#define DEFAULT_ADS1115_CONFIGS 0b0000'0100'1000'0011

/** @brief Gère l'interfaçage avec l'ADC externe en I2C
 */
class ADS1115 : private I2C {
public:
  uint16_t configs;

  ADS1115(const uint8_t deviceAddr);

  /** @brief Assigne le bon channel sur l'ADC
   *  @param channel Le channel à lire entre 0 et 3
   */
  void setChannel(const uint8_t channel);

  /** @brief Lis une valeur sur un channel de l'ADC
   *  @return La valeur 16 bits lu
   */
  uint16_t read();
};

#endif // ADS_1115_H
