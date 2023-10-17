#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath,
         std::array<uint8_t, 4> pins)
    : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins) {
  state = Initialize;
}

void STM::start(){
  switch (state){
    case Initialize:
      initialize();
      // state = Find_sample;
      break;
    case Find_sample:

      break;
    case Lower_motor:
      break;
    case Mesure_height:
      break;
    case Save_pixel:
      break;
    case Goto_next_coordinate:
      break;
  }
}

void STM::initialize() {
  StepMotor::goToRelative(-40);
}

void STM::acquirePixelAtConstantHeight() {}

void STM::acquirePixelAtConstantCurrent() {}

void STM::goToNextPosition() {}

void STM::fixImage() {}

void STM::exportImage(const std::string filename) {}
