#include <cstdint>
#include <iostream>

#include "I2C.h"
#include "Plateforme.h"
#include "Utils.h"
#include "Vector3D.h"

int main() {
  Plateforme plateforme("/dev/ttyUSB0");
  I2C conn(10);
}
