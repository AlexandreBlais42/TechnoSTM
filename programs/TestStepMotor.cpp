#include <iostream>
#include <array>
#include <cstdint>

#include "StepMotor.h"

int main (int argc, char *argv[]) {
std::array<uint8_t, 4> pins = {6, 7, 8, 9};
  StepMotor motor(0, pins);
  for (const auto pin : motor.pins){
    std::cout << std::to_string(pin) << std::endl;
  }
  motor.goToRelative(-5);
  return 0;
}
