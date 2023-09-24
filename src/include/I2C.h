#ifndef I2C_H
#define I2C_H

#include <cstdint>
#include <iostream>

#ifdef __arm__
#include <wiringPiI2C.h>
#endif // __arm__

/** @brief Classe qui s'occupe de la communication I2C
 */
class I2C {
public:
  int fileDescriptor;

  I2C(const uint8_t deviceAddr);
  ~I2C();

  /** @brief Écrit une valeur de 8 bit sur un registre
   *  @param reg L'addresse du registre à écrire
   *  @param data La donnée 8 bits à écrire
   */
  void writeReg8(int reg, uint8_t data);
  /** @brief Écrit une valeur de 16 bit sur un registre
   *  @param reg L'addresse du registre à écrire
   *  @param data La donnée 16 bits à écrire
   */
  void writeReg16(int reg, uint16_t data);
  /** @brief Lit une valeur 8 bits d'un registre
   *   @param reg L'addresse du registre à lire
   *   @return La valeur 8 bit lu
   */
  uint8_t readReg8(int reg);
  /** @brief Lit une valeur 16 bits d'un registre
   *   @param reg L'addresse du registre à lire
   *   @return La valeur 16 bit lu
   */
  uint16_t readReg16(int reg);
};

#endif // I2C_H
