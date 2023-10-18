#include "Aiguille.h"

Aiguille::Aiguille(const uint8_t deviceAddr) : ADS1115(deviceAddr) {}

int16_t Aiguille::readVoltage() { return read(); }
