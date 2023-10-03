#include "STM.h"

STM::STM(const uint8_t deviceAddr, const std::string devicePath, std::array<uint8_t, 4> pins) : Aiguille(deviceAddr), Plateforme(devicePath), StepMotor(0, pins){

}
