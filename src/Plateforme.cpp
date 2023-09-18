#include "Plateforme.h"

Plateforme::Plateforme(const std::string devicePath){
  ser = Serial();
  ser.begin(devicePath, 115200);
}

inline void Plateforme::goToRelative(const uint16_t x, const uint16_t y, const uint16_t z){
  position += Vector3D<uint16_t>(x, y, z);
}

inline void Plateforme::goToRelative(const Vector3D<uint16_t> &coordinates){
  position += coordinates;
}


inline void Plateforme::goToAbsolute(const uint16_t x, const uint16_t y, const uint16_t z){
  position = Vector3D<uint16_t>(x, y, z);
}

inline void Plateforme::goToAbsolute(const Vector3D<uint16_t> &coordinates){
  position = coordinates;
}

void Plateforme::goToPosition() const {
  std::vector<std::string> stringsToSend(0);

  ser.write("A" + getHexString(position.x));
  ser.write("B" + getHexString(position.y));
  ser.write("C" + getHexString(position.z));
}
