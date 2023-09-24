#include "Serial.h"

Serial::Serial() {}

void Serial::begin(const std::string devicePath, const uint32_t baudrate) {
#ifdef __arm__
  // À implémenter
#else
  std::cout << "Serial::begin, device : " << devicePath
            << " baudrate : " << baudrate << std::endl;
#endif // __arm__
}

void Serial::write(const std::string s) const {
#ifdef __arm__
  // À implémenter
#else
  std::cout << "Serial::write : " << s << std::endl;
#endif // __arm__
}

std::vector<char> Serial::read(const uint32_t amount) const {
#ifdef __arm__
  // À implémenter
  return std::vector<char>(0);
#else
  std::cout << "Serial::read : " << amount << std::endl;
  return std::vector<char>(0);
#endif // __arm__
}
