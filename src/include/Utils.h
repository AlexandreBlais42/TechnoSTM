#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

/** @brief Retourne l'équivalent hexadécimal de l'argument
 *  @param val La valeur plus petite que 16 à convertir
 *  @return Le charactère converti
 */
char getHexCharacter(const uint8_t val);

/** @brief Obtiens la représentation hexadécimal du de la donnée passée
 *  @param data La donnée à convertir
 *  @return Chaîne de charactères hexadécimal
 */
std::string getHexString(uint16_t data);

void delay_ms(const uint32_t time);

/** @brief Inverse l'ordre des octets dans un type de donnée quelconque
 *  @param data Les données à inverser
 *  @return Les données inversées
 *  @example 0x12345678 deviendrait 0x78563412
 */
template <typename T> T invertBytes(T data) {
  T toReturn = 0;
  uint8_t *ptr = (uint8_t *)&data;
  for (uint16_t _ = 0; _ < sizeof(data); _++) {
    toReturn <<= 8;
    toReturn += *ptr;
    ptr++;
  }
  return toReturn;
}

template <typename T> std::vector<uint8_t> getBytes(T data) {
  uint8_t *bytes = reinterpret_cast<uint8_t *>(&data); 
  std::vector<uint8_t> toReturn(0);
  toReturn.reserve(sizeof(data));
  toReturn.insert(toReturn.begin(), bytes, bytes + sizeof(data));
  return toReturn;
}

void writeStringToVector(const std::string s, std::vector<uint8_t> v);

#endif // UTILS_H
