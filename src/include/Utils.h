#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <chrono>
#include <iostream>
#include <string>
#include <thread>

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
 *  @example 0x1234 deviendrait 0x3412
 */
template<typename T>
T invertBytes(T data);

#endif // UTILS_H
