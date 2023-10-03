#ifndef STM_H
#define STM_H

#include <iostream>
#include <cstdint>

#include "Aiguille.h"
#include "Plateforme.h"
#include "StepMotor.h"

/** @brief Classe qui gère les fonctions du microscope
 */
class STM : private Aiguille, private Plateforme, private StepMotor {
public:
  STM(const uint8_t deviceAddr, const std::string devicePath, std::array<uint8_t, 4> pins);


};

#endif // STM_H
