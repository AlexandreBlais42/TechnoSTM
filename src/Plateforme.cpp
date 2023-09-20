#include "Plateforme.h"

Plateforme::Plateforme(const std::string devicePath){
  ser = Serial();
  ser.begin(devicePath, 115200);
}

inline void Plateforme::setPositionRelative(const uint16_t x, const uint16_t y, const uint16_t z){
  position += Vector3D<uint16_t>(x, y, z);
}

inline void Plateforme::setPositionRelative(const Vector3D<uint16_t> &coordinates){
  position += coordinates;
}


inline void Plateforme::setPositionAbsolute(const uint16_t x, const uint16_t y, const uint16_t z){
  position = Vector3D<uint16_t>(x, y, z);
}

inline void Plateforme::setPositionAbsolute(const Vector3D<uint16_t> &coordinates){
  position = coordinates;
}

void Plateforme::moveToPosition() const {
  std::vector<std::string> stringsToSend(0);

  ser.write("A" + getHexString(position.x));
  ser.write("B" + getHexString(position.y));
  ser.write("C" + getHexString(position.z));
}
