#include <cstdint>
#include <iostream>

#include "STM.h"

int main(int argc, char *argv[]) {
  std::array<uint8_t, 4> pinsStepMotor = {25, 8, 7, 1};
  STM microscope(0x48, "/dev/ttyUSB0", pinsStepMotor);
  microscope.start();
}
