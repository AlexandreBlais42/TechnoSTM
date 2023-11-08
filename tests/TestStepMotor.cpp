#include <array>
#include <cstdint>
#include <iostream>

#include "StepMotor.h"

int main(int argc, char *argv[]) {
  int32_t numberOfSteps = -500; // Valeur par d�faut
  if (argc >= 2){
    numberOfSteps = std::stoi(argv[1]);
  }
  std::array<uint8_t, 4> pins = {25, 8, 7, 1};
  StepMotor motor(0, pins);
  for (const auto pin : motor.pins) {
    std::cout << std::to_string(pin) << std::endl;
  }
  motor.setPositionRelative(numberOfSteps);
  return 0;
}
