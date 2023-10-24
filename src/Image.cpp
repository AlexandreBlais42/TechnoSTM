#include "Image.h"

Image::Image(uint32_t XRes, uint32_t YRes, std::string title, std::string ZUnits) : XRes(XRes), YRes(YRes), title(title), ZUnits(ZUnits){
  XReal = XRes * 0.000'000'000'15;
  YReal = YRes * 0.000'000'000'15;
  XOffset = 0;
  YOffset = 0;
  XYUnits = "m";
  data = std::vector<std::vector<float>>(XRes, std::vector<float>(YRes)); // Empty Vector
}

std::vector<uint8_t> Image::getGSFFile(){
  std::vector<uint8_t> GSFFile;
  GSFFile.reserve(XRes * YRes * 4 + 1024); // Réserver assez d'espace pour stocker l'entièreté du fichire
  std::string header = 
    std::string("Gwyddion Simple Field 1.0") + 
    "\nXRes = " + std::to_string(XRes) + 
    "\nYRes = " + std::to_string(YRes) + 
    "\nXReal = 15e-8" + 
    "\nYReal = 15e-8" + 
    "\nXYUnits = " + XYUnits +
    "\nZUnits = " + ZUnits + 
    "\nTitle = " + title;
    
  writeStringToVector(header, GSFFile);
  for (uint8_t i = 0 ; i < 4 - (header.size() % 4) ; i++){
    GSFFile.push_back('\0');
  }

  for (std::vector<float> row : data){
    for (const float &value : row){
      for (const uint8_t &byte : getBytes(value)){
        GSFFile.push_back(byte);
      }
    }
  }

  return GSFFile;
}
