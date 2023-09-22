#include <iostream>
#include <cstdint>

#include "Plateforme.h"
#include "Vector3D.h"
#include "Utils.h"
#include "I2C.h"

int main(){
  Plateforme plateforme("/dev/ttyUSB0");
  I2C conn(10);
}
