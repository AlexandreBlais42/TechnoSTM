#include <array>
#include <cstdint>
#include <iostream>

#include "StepMotor.h"

int main(int argc, char *argv[]) {
  std::array<uint8_t, 4> pins = {25, 8, 7, 1};
  StepMotor motor(0, pins);
  for (const auto pin : motor.pins) {
    std::cout << std::to_string(pin) << std::endl;
  }
  motor.goToRelative(-500);
  return 0;
}
