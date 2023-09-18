#include "Serial.h"

Serial::Serial(){}

void Serial::begin(const std::string devicePath, const uint32_t baudrate){
  // À implémenter
  std::cout << "Serial::begin, device : " << devicePath << " baudrate : " << baudrate << std::endl;
}

void Serial::write(const std::string s) const {
  // À implémenter
  std::cout << "Serial::write : " << s << std::endl;
}

std::vector<char> Serial::read(const uint32_t amount) const {
  //À implémenter
  std::cout << "Serial::read : " << amount << std::endl;
  return std::vector<char>(0);
}
