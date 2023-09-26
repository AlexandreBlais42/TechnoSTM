#include "I2C.h"
#include <iostream>
/*  Ce fichier ce décompose en deux parties: le code pour le raspberry pi
 * (#ifdef __arm__) et le code pour les autres ordi. Le code doit être séparé
 * car les ordinateurs n'ont pas accès au GPIO comme le raspberry pi et il
 * serait impossible de le compiler sinon.
 */
#ifdef __arm__

I2C::I2C(const uint8_t deviceAddr) {
  fileDescriptor = wiringPiI2CSetup(deviceAddr);
  if (!fileDescriptor) {
    std::cout << "wiringPiI2CSetup n'a pas pu être initialisé" << std::endl;
    exit(1);
  }
}

I2C::~I2C() {
  // À implémenter
}

void I2C::writeReg8(int reg, uint8_t data) {
  wiringPiI2CWriteReg8(fileDescriptor, reg, data);
}
void I2C::writeReg16(int reg, uint16_t data) {
  wiringPiI2CWriteReg16(fileDescriptor, reg, data);
}
uint8_t I2C::readReg8(int reg) {
  return wiringPiI2CReadReg8(fileDescriptor, reg);
}
uint16_t I2C::readReg16(int reg) {
  return wiringPiI2CReadReg16(fileDescriptor, reg);
}

#else // Pas sur le raspberry pi

I2C::I2C(const uint8_t deviceAddr) {
  std::cout << "I2C::I2C appellé avec deviceAddr = "
            << std::to_string(deviceAddr) << std::endl;
}

I2C::~I2C() { std::cout << "I2C::~I2C appellé " << std::endl; }

void I2C::writeReg8(int reg, uint8_t data) {
  std::cout << "I2C::" << __FUNCTION__ << "appellé avec reg = " << reg
            << " et data = " << data << std::endl;
}
void I2C::writeReg16(int reg, uint16_t data) {
  std::cout << "I2C::" << __FUNCTION__ << "appellé avec reg = " << reg
            << " et data = " << data << std::endl;
}

uint8_t I2C::readReg8(int reg) {
  std::cout << "I2C::" << __FUNCTION__ << "appellé avec reg = " << reg
            << std::endl;
  return 0;
}

uint16_t I2C::readReg16(int reg) {
  std::cout << "I2C::" << __FUNCTION__ << "appellé avec reg = " << reg
            << std::endl;
  return 0;
}

#endif // __arm__
