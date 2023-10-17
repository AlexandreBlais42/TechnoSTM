#include <cstdint>
#include <iostream>

#include "STM.h"

int main(int argc, char *argv[]) {
  std::array<uint8_t, 4> pins = {25, 8, 7, 3};
  STM microscope(0x48, "/dev/ttyUSB0", pins);
}
