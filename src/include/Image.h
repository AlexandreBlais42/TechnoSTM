#ifndef IMAGE_H
#define IMAGE_H

#include <cstdint>
#include <iostream>
#include <string> 
#include <vector>

#include "Utils.h"

/** @brief Classe pour l'image et le traitement de celle ci et s'exporte en format gsf (https://gwyddion.net/documentation/user-guide-en/gsf.html)
 */
class Image {
public:
  // Ces variables viennent directement du format gsf de Gwyddion
  uint32_t XRes, YRes;
  float XReal, YReal;
  float XOffset, YOffset;
  std::string title;
  std::string XYUnits;
  std::string ZUnits;
  std::vector<std::vector<float>> data;

  Image(const uint32_t XRes, const uint32_t YRes, const std::string title, const std::string ZUnits);

  /** @brief Cette fonction retourne un fichier en format *.gsf dans un vecteur de uint8_t
   *  @return Le vecteur d'uint8_t avec les données
   */
  std::vector<uint8_t> getGSFFile();
};

#endif // IMAGE_H
