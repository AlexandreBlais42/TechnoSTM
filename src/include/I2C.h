#ifndef I2C_H
#define I2C_H

#include <cstdint>
#include <iostream>
#include <unistd.h>

#include "Utils.h"
#include "Serial.h"

#ifdef __arm__
/** @note À changer pour <linux-i2c-dev.h> si ça marche pas bien
 */
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

  template <typename T>
  /** @brief Écrits une donnée sur la communication I2C
   */
  void writeData(const T &data) {
#ifdef __arm__
    T invertedData = invertBytes<T>(data);
    if (!write(fileDescriptor, (char *)&invertedData, sizeof(T))) {
      std::cout << "Une erreur est survenue lors de l'écriture dans la "
                   "communication I2C";
    }
#else
    DEBUG << "I2C::writeData appellé avec data (" << sizeof(data) << " octets)"
              << std::endl;
#endif // __arm__
  }

  template <typename T>
  /** @brief Lits une donnée sur la communication I2C
   */
  T readData() {
#ifdef __arm__
    T data;
    read(fileDescriptor, (char *)&data, sizeof(T));
    data = invertBytes<T>(data);
    return data;
#else
    DEBUG << "I2C::readData appellé" << std::endl;
    return 0;
#endif // __arm__
  }
};

#endif // I2C_H
