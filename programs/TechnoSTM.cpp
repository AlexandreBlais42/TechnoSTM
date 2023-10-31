#include <cstdint>
#include <iostream>
#include <queue>

#include "STM.h"
#include "Error.h"

int main(int argc, char *argv[]) {
  std::array<uint8_t, 4> pinsStepMotor = {25, 8, 7, 1};
  STM microscope(0x48, "/dev/ttyUSB0", pinsStepMotor, 500, 500, 1);
  std::queue<uint8_t> controlData;
  microscope.start();
}
