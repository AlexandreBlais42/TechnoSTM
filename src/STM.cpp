#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath,
         std::array<uint8_t, 4> pins)
    : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins) {}

void STM::initialize() {}

void STM::acquirePixelAtConstantHeight() {}

void STM::acquirePixelAtConstantCurrent() {}

void STM::goToNextPosition() {}

void STM::fixImage() {}

void STM::exportImage(const std::string filename) {}
