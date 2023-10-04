#include "Serial.h"

Serial::Serial(const std::string devicePath, const uint32_t baudrate) {
  begin(devicePath, baudrate);
}

void Serial::begin(const std::string devicePath, const uint32_t baudrate) {
#ifdef __arm__
  if (fileDescriptor = serialOpen(devicePath.c_str(), baudrate) == -1) {
    std::cout << "Erreur lors de l'initialisation de la connection sérielle"
              << std::endl;
    exit(1);
  }
#else
  std::cout << "Serial::begin, device : " << devicePath
            << " baudrate : " << baudrate << std::endl;
#endif // __arm__
}

void Serial::write(const std::string s) const {
#ifdef __arm__
  serialPuts(fileDescriptor, s.c_str());
#else
  std::cout << "Serial::write : " << s << std::endl;
#endif // __arm__
}

std::vector<char> Serial::read(const uint32_t amount) const {
#ifdef __arm__
  std::vector<char> toReturn(0);
  for (uint32_t i = 0; i < amount; i++) {
    toReturn.push_back(serialGetchar(fileDescriptor));
  }
#else
  std::cout << "Serial::read : " << amount << std::endl;
  return std::vector<char>(amount, 0);
#endif // __arm__
}
