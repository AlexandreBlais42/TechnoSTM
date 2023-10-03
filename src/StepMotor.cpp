#include "StepMotor.h"

StepMotor::StepMotor(const int32_t position, const std::array<uint8_t, 4> pins) : position(position), pins(pins), pinIndex(0) {
  for (const uint8_t pin : pins){
    pinMode(pin, OUTPUT);
  }
}

void StepMotor::goToRelative(const int32_t pos) {
  position += pos;

  int8_t direction = 1;
  if (pos < 0) {
    direction *= -1;
  }

  for (uint8_t _ = 0; _ < abs(pos); _++) {
    setPin(pins[pinIndex], false);
    pinIndex += direction;
    pinIndex += 4; // -1 % 4 == -1 en C++, On a besoin d'un index positif
    pinIndex %= 4;
    setPin(pins[pinIndex], true);

    delay_ms(5);
  }
}

inline void StepMotor::goToAbsolute(const int32_t pos) {
  goToRelative(pos - position);
}
