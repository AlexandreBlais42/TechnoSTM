#include <iostream>
#include <cstdint>

#include "Plateforme.h"
#include "Vector3D.h"

int main(){
  Plateforme plateforme("/dev/ttyUSB0");

  plateforme.goToPosition();
}
