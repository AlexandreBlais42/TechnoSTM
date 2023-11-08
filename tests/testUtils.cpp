#include <cassert>
#include <iostream>
#include <vector>

#include "Utils.h"

int main(int argc, char *argv[]) {

  // Test invertBytes
  std::cout << "Test invertBytes" << std::endl;
  assert(invertBytes<uint32_t>(0x12345678) == 0x78563412);

  // Test getBytes
  std::cout << "Test getBytes" << std::endl;
  std::vector<uint8_t> bytes = getBytes<uint32_t>(0x00000101);
  assert(bytes[0] == 1);
  assert(bytes[1] == 1);
  assert(bytes[2] == 0);
  assert(bytes[3] == 0);

  std::cout << "Tous les tests ont étés passés" << std::endl;
  return 0;
}
