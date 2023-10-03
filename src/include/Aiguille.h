#ifndef AIGUILLE_H
#define AIGUILLE_H

#include <cstdint>

#include "ADS1115.h"

/** @brief Classe qui gère l'aiguille
 */
class Aiguille : private ADS1115 {
public:
  /** @brief Constructeur de la classe aiguille
  *   @param ADCAddr
  */
  Aiguille(const uint8_t deviceAddr);

  /** @brief Lit une valeur correspondant au courant sur l'aiguille
  *   @return La valeur du courant lu
  */
  int16_t readVoltage();
};

#endif // AIGUILLE_H
