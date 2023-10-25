#include "StepMotor.h"

StepMotor::StepMotor(const int32_t position, const std::array<uint8_t, 4> pins)
    : position(position), pins(pins), pinIndex(0) {
  for (const uint8_t pin : pins) {
    _pinMode(pin, OUTPUT);
  }
}

void StepMotor::setPositionRelative(const int32_t pos) {
  position += pos;

  int8_t direction = 1;
  if (pos < 0) {
    direction *= -1;
  }

  for (uint16_t _ = 0; _ < abs(pos); _++) {
    setPin(pins[pinIndex], false);
    pinIndex += direction;
    pinIndex %= 4;
    setPin(pins[pinIndex], true);

    delay_ms(3); // @note Le Délai optimal est à trouver
  }
}

inline void StepMotor::setPositionAbsolute(const int32_t pos) {
  setPositionRelative(pos - position);
}
