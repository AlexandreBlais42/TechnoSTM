#ifndef ADS1115_H
#define ADS1115_H

#include <cassert>
#include <cstdint>
#include <iostream>

#include "I2C.h"

#define DEFAULT_ADS1115_CONFIGS 0b0000'0100'1000'0011

#define OPERATIONAL_STATUS_MASK 0b0111'1111'1111'1111
#define OPERATIONAL_STATUS_OFFSET 15
#define CHANNEL_MASK 0b1000'1111'1111'1111
#define CHANNEL_OFFSET 12
#define PROGRAMMABLE_GAIN_MASK 0b1111'0001'1111'1111
#define PROGRAMMABLE_GAIN_OFFSET 9
#define SINGLE_CONVERSION_MASK 0b1111'1110'1111'1111
#define SINGLE_CONVERSION_OFFSET 8
#define DATA_RATE_MASK 0b1111'1111'0001'1111
#define DATA_RATE_OFFSET 5
#define COMPARE_MODE_MASK 0b1111'1111'1110'1111
#define COMPARE_MODE_OFFSET 4
#define COMPARATOR_POLARITY_MASK 0b1111'1111'1111'0111
#define COMPARATOR_POLARITY_OFFSET 3
#define LATCHING_COMPARATOR_MASK 0b1111'1111'1111'1011
#define LATCHING_COMPARATOR_OFFSET 2
#define COMPARATOR_QUEUE_MASK 0b1111'1111'1111'1100
#define COMPARATOR_QUEUE_OFFSET 0

/** @brief Gère l'interfaçage avec l'ADC externe en I2C
 *  @note https://www.ti.com/lit/ds/symlink/ads1115.pdf
 */
class ADS1115 : private I2C {
private:
  /** @brief gère l'addressPointer
   */
  typedef enum AddressPointer_t {
    Conversion = 0b00,
    Config = 0b01,
    Lo_thresh = 0b10,
    Hi_thresh = 0b11,
  } AddressPointer_t;

public:
  uint16_t configs;
  AddressPointer_t addressPointer; // Page 27 dans la documentation

  ADS1115(const uint8_t deviceAddr);

  /** @brief Assigne le bit operational status
   *  @param status Vrai: start conversion Faux: no effect
   */
  void setOperationalStatus(const bool status);

  /** @brief Assigne le bon channel sur l'ADC
   *  @param channel Deux bits qui représentent les channels de 0 à 3 inclusif
   */
  void setChannel(const uint8_t channel);

  /** @brief Assigne le gain conformément à la datasheet
   *  @param gain Une valeur de 3bits défini à la page 28 de la datasheet
   *  @note Dans la datasheet, ça affecte le FSR
   */
  void setProgrammableGain(const uint8_t gain);

  /** @brief Assigne le bit de single conversion
   *  @param mode Vrai: Single-conversion, Faux: Continuous-conversion
   */
  void setSingleConversionMode(const bool mode);

  /** @brief Assigne le datarate
   *  @param rate Une valeur de 3 bits spécifiant le data rate
   */
  void setDataRate(const uint8_t rate);

  /** @brief Assigne le mode de comparaison
   *  @param mode Vrai: window comparator, Faux: traditional comparator
   */
  void setCompareMode(const bool mode);

  /** @brief Assigne la polarité du comparateur
   *  @param mode Vrai: Active high, Faux: Active low
   */
  void setComparatorPolarity(const bool mode);

  /** @brief Assigne le latch de la pin ALERT/RDY
   *  @param mode Vrai: La pin reste latch jusqu'à la lecture de la donnée,
   * Faux: Pas de latch
   */
  void setLatchingComparator(const bool mode);

  /** @brief Assigne le mode de fifo du comparateur
   *  @param mode Le mode à mettre, voir la datasheet
   */
  void setComparatorQueue(const uint8_t mode);

  /** @brief Écris les configs sur l'ADS1115
   */
  void writeConfigs();

  /** @brief Lis une valeur sur un channel de l'ADC
   *  @return La valeur 16 bits lu
   */
  int16_t read();
};

#endif // ADS_1115_H
