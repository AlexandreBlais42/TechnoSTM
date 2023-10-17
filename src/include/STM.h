#ifndef STM_H
#define STM_H

#include <cstdint>
#include <iostream>

#include "Aiguille.h"
#include "Plateforme.h"
#include "StepMotor.h"

// @todo Trouver la valeur de lecture d'aiguille à laquelle arrêter
#define AIGUILLE_THRESHOLD_VOLTAGE 1000             // Voltage à lequel on considère que le matériel a été détecté
#define AIGUILLE_CONSTANT_CURRENT_VOLTAGE 20000     // Voltage (arbitraire) que l'on doit avoir en mode courant constant


/** @brief Classe qui gère les fonctions du microscope
 */
class STM : private Aiguille, private Plateforme, private StepMotor {
public:
  /** @brief États de la machine à états du microscope
   *  @note Initialize            : Monte le step moteur de 40 pour s'éloigner du matériel
   *  @note Find_sample           : Monte le matériel jusqu'à ce qu'il soit détecté ou que la hauteur maximale est atteinte
   *  @note Lower_motor           : Descend le matériel et ensuite le moteur
   *  @note Mesure_height         : Ajuste le microscope jusqu'à obtenir le voltage désiré
   *  @note Save_pixel            : Insère la hauteur mesurée dans l'image
   *  @note Goto_next_coordinate  : Déplace le microscope à la prochaine coordonnée pour l'image
   */
  typedef enum StateMachineStates {
    Initialize,
    Find_sample,
    Lower_motor,
    Mesure_height,
    Save_pixel,
    Goto_next_coordinate
  } StateMachineStates;

  STM(const uint8_t deviceAddr, const std::string devicePath,
      std::array<uint8_t, 4> pins);

  StateMachineStates state;

  /** @brief Démarre la machine à état qui gère le microscope
   */
  void start();

private:
  /** @brief Lit le voltage de l'aiguille et l'insère dans l'image (Mesure par
   * hauteur constante)
   */
  void acquirePixelAtConstantHeight();

  /** @brief Lit la hauteur de l'aiguille pour un certain courant et l'insère
   * dans l'image (Mesure par courant constant)
   */
  void acquirePixelAtConstantCurrent();

  /** @brief Va à la prochaine coordonnée pour la prochaine mesure
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
