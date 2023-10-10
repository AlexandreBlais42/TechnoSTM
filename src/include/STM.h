#ifndef STM_H
#define STM_H

#include <cstdint>
#include <iostream>

#include "Aiguille.h"
#include "Plateforme.h"
#include "StepMotor.h"

/** @brief Classe qui gère les fonctions du microscope
 */
class STM : private Aiguille, private Plateforme, private StepMotor {
public:
  STM(const uint8_t deviceAddr, const std::string devicePath,
      std::array<uint8_t, 4> pins);

  /** @brief Initialise la position de la plateforme et du stepmoteur
   */
  void initialize();

  /** @brief Lit le voltage de l'aiguille et l'insère dans l'image (Mesure par
   * hauteur constante)
   */
  void acquirePixelAtConstantHeight();

  /** @brief Lit la hauteur de l'aiguille pour un certain courant et l'insère
   * dans l'image (Mesure par courant constant)
   */
  void acquirePixelAtConstantCurrent();

  /** @brief Va à la prochaine coordonée pour la prochaine mesure
   */
  void goToNextPosition();

  /** @brief Enlève la composante exponentielle de l'image dû au fait que le
   * courant de tunneling dépend de la distance de manière non linéaire? (au
   * carré ou exponentielle ?)
   * @todo Trouver fonction qui décrit le courant de tunneling en fonction de la
   * distance
   * @note On va surement utiliser Gwyddion pour le traitement d'image, à voir.
   */
  void fixImage();

  /** @brief Sauvegarde l'image dans un fichier
   */
  void exportImage(const std::string filename);
};

#endif // STM_H
