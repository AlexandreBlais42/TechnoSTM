#ifndef SERIAL_H
#define SERIAL_H

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "Vector3D.h"

/** @brief Classe qui permet de communiquer sériellement entre autre pour
 * communiquer avec le DAC 16bits AnalogDevices
 */
class Serial {
public:
  std::string devicePath;

  Serial();

  /** @brief Setup la connection sérielle
   *  @param devicePath Le chemin du fichier qui permet de communiquer sur le
   * port sériel
   *  @param baudrate Le débit en bit/s de la communication
   */
  void begin(const std::string devicePath, const uint32_t baudrate);

  /** @brief Écrit sur le port sériel
   *  @param s La chaîne de charactères à écrire
   */
  void write(const std::string s) const;

  /** @brief Lis des données du port sériel
   *  @param amount La quantité d'octets à lire
   *  @return Un vecteur de char contenant les données lues
   */
  std::vector<char> read(const uint32_t amount) const;
};

#endif // SERIAL_H
