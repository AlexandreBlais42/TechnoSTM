#include <cstdint>
#include <cstring>
#include <iostream>

#include "Serial.h"
#include "Utils.h"

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << "Usage : " << argv[0]
              << "/dev/ttyUSBx channel voltage(en hexadécimal)" << std::endl;
    return EXIT_FAILURE;
  }

  Serial connection(argv[1], 115200);
  connection.write(std::string(argv[2]) + std::string(argv[3]));
  return 0;
}
