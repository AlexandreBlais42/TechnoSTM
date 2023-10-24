#include "Image.h"

Image::Image(){}

void Image::initialize(const uint32_t XRes, const uint32_t YRes, const std::string title, const std::string ZUnits){
  this->XRes = XRes;
  this->YRes = YRes;
  this->title = title;
  this->ZUnits = ZUnits;
  XReal = XRes * 0.000'000'000'15;
  YReal = YRes * 0.000'000'000'15;
  XOffset = 0;
  YOffset = 0;
  XYUnits = "m";
  data = std::vector<std::vector<float>>(XRes, std::vector<float>(YRes)); // Empty Vector
}

void Image::setPixel(const uint32_t XCoordinate, const uint32_t YCoordinate, const float value){
  assert(XCoordinate < XRes);
  assert(YCoordinate < YRes);
  data[YCoordinate][XCoordinate] = value;
}

std::vector<uint8_t> Image::getGSFFile(){
  // https://gwyddion.net/documentation/user-guide-en/gsf.html
  std::vector<uint8_t> GSFFile;

  // Réserver assez d'espace dans le vecteur pour ne pas devoir en réserver dans chaque push_back
  GSFFile.reserve(XRes * YRes * 4 + 1024);

  // Écrits le header bien formatté
  std::string header = 
    std::string("Gwyddion Simple Field 1.0") + 
    "\nXRes = " + std::to_string(XRes) + 
    "\nYRes = " + std::to_string(YRes) + 
    "\nXReal = 15e-11" + 
    "\nYReal = 15e-11" + 
    "\nXYUnits = " + XYUnits +
    "\nZUnits = " + ZUnits + 
    "\nTitle = " + title;
    
  // Ajoute les NULL de pad
  writeStringToVector(header, GSFFile);
  for (uint8_t i = 0 ; i < 4 - (header.size() % 4) ; i++){
    GSFFile.push_back('\0');
  }

  // Rajoute chaque données dans le vecteur
  for (std::vector<float> row : data){
    for (const float &value : row){
      for (const uint8_t &byte : getBytes(value)){
        GSFFile.push_back(byte);
      }
    }
  }

  GSFFile.shrink_to_fit();

  return GSFFile;
}
