#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <cassert>

/** @brief Retourne l'équivalent hexadécimal de l'argument
 *  @param val La valeur plus petite que 16 à convertir
 *  @return Le charactère converti
 */
char map4BitToHexChar(const uint8_t val);

/** @brief Obtiens la représentation hexadécimal du de la donnée passée
 *  @param data La donnée à convertir
 *  @return Chaîne de charactères hexadécimal
 */
std::string getHexString(uint16_t data);

void delay_ms(const uint32_t time);

#endif //UTILS_H
