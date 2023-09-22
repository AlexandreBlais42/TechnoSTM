#ifndef ADS1115_H
#define ADS1115_H

#include <iostream>
#include <cstdint>
#include <cassert>
/** @brief Gère l'interfaçage avec l'ADC externe en I2C
 */
class ADS1115{
public:
  ADS1115();

  /** @brief Lis une valeur sur un channel de l'ADC
   *  @param channel Le channel à lire entre 0 et 3
   *  @return La valeur 16 bits lu
   */
  uint16_t read(const uint8_t channel) const ;
  
};

#endif // ADS_1115_H
