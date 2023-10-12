#include <cstdint>
#include <iostream>

#include "Serial.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage : " << argv[0] << "/dev/ttyUSBx" << std::endl;
    return EXIT_FAILURE;
  }

  Serial connection(argv[1], 115200);
  connection.write("Test d'ectriture");

  std::cout << "Lecture du port sériel : ";
  for (char c : connection.read(5)) {
    std::cout << c;
  }
  std::cout << std::endl;
  return 0;
}
