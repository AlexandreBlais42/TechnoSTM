#include "Plateforme.h"

Plateforme::Plateforme(const std::string devicePath)
    : Serial(devicePath, 115200) {
  Serial::write("D00FF\n"); // Environ 50mV pour le voltage de biais
}

void Plateforme::setPositionRelative(const uint16_t x, const uint16_t y,
                                     const uint16_t z) {
  position += Vector3D<uint16_t>(x, y, z);
  moveToPosition();
}

void Plateforme::setPositionRelative(const Vector3D<uint16_t> &coordinates) {
  position += coordinates;
  moveToPosition();
}

void Plateforme::setPositionAbsolute(const uint16_t x, const uint16_t y,
                                     const uint16_t z) {
  position = Vector3D<uint16_t>(x, y, z);
  moveToPosition();
}

void Plateforme::setPositionAbsolute(const Vector3D<uint16_t> &coordinates) {
  position = coordinates;
  moveToPosition();
}

void Plateforme::moveToPosition() const {
  Serial::write("A" + getHexString(position.x) + "\n");
  Serial::write("B" + getHexString(position.y) + "\n");
  Serial::write("C" + getHexString(position.z) + "\n");
}
