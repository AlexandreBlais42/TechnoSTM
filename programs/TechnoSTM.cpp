#include <cstdint>
#include <iostream>

#include "STM.h"
#include "Debug.h"

int main(int argc, char *argv[]) {
  DEBUG << "test" << std::endl;
  std::array<uint8_t, 4> pinsStepMotor = {25, 8, 7, 1};
  STM microscope(0x48, "/dev/ttyUSB0", pinsStepMotor);
  microscope.start();
}
