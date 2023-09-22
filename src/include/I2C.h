#ifndef I2C_H
#define I2C_H

#include <iostream>
#include <fstream>
#include <cstdint>

#include <linux/i2c-dev.h>
/*  @todo Utiliser <wiringPi.h>
 */

/** @brief Classe qui s'occupe de la communication I2C
 */
class I2C{
public:
  std::fstream file; 

  I2C();
  I2C(const uint8_t deviceAddr);
  ~I2C();

  template<typename T>
  /** @brief Écrits une donnée sur la communication I2C
   */
  void write(const T &data){
    file.write((char *) &data, sizeof(T));
  }

  template<typename T>
  /** @brief Lits une donnée sur la communication I2C
   */
  T read(){
    T data;
    file.read((char *) &data, sizeof(T));
    return data;
  }
};

#endif // I2C_H
