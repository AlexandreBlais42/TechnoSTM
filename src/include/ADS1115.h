#ifndef ADS1115_H
#define ADS1115_H

#include <cassert>
#include <cstdint>
#include <iostream>

#include "I2C.h"

#define DEFAULT_ADS1115_CONFIGS 0b0000'0100'1000'0011
#define CHANNEL_MASK 0b1000'1111'1111'1111

/** @brief Gère l'interfaçage avec l'ADC externe en I2C
 *  @note https://www.ti.com/lit/ds/symlink/ads1115.pdf
 */
class ADS1115 : private I2C {
private:
  /** @brief gère l'addressPointer
   */
  typedef enum AddressPointer_t {
    Conversion,
    Config,
    Lo_thresh,
    Hi_thresh,
  } AddressPointer_t;

public:
  uint16_t configs;
  AddressPointer_t addressPointer; // Page 27 dans la documentation

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
