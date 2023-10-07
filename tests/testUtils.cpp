#include <cassert>
#include <iostream>

#include "Utils.h"

int main(int argc, char *argv[]) {

  uint32_t a = invertBytes<uint32_t>(0x12345678);
  std::cout << a << std::endl;
  return 0;
}
