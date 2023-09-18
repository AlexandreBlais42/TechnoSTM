#include "StepMotor.h"

StepMotor::StepMotor(const int32_t position) : position(position), steps(0){

}


void StepMotor::goToRelative(const int32_t pos){
  position += pos;

  int8_t direction = 1;
  if (pos < 0){
    direction *= -1;
  }

  for (uint8_t _ = 0 ; _ < abs(pos) ; _++){
    GPIO::setPin(steps, false);
    steps += direction;
    steps %= 4;
    GPIO::setPin(steps, true);

    delay_ms(5);
  }
}

inline void StepMotor::goToAbsolute(const int32_t pos){
  goToRelative(pos - position);
}
