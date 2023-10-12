#include <iostream>

#include "ADS1115.h"

int main(int argc, char *argv[]) {
  std::cout << "Test de l'ADC, appuyer sur enter pour lire une valeur";
  ADS1115 adc(0x48);
  adc.setChannel(0);
  adc.setSingleConversionMode(false);
  adc.writeConfigs();

  int16_t readVoltage;
  while (true) {
    std::cin.get(); // Attends un enter pour lire le voltage
    readVoltage = adc.read();

    // Le 2.048 vient du FSR par défaut
    std::cout << "Valeur lu : " << readVoltage << " ("
              << (float)readVoltage / 32768 * 2.048 << " V)";
  }
}
