#include "Utils.h"

char getHexCharacter(const uint8_t val) {
  assert(val < (1 << 4));
  if (val <= 9) {
    return '0' + val;
  }
  return 'A' + val - 10;
}

std::string getHexString(uint16_t data) {
  std::string str = "";
  char chars[3];
  for (uint8_t i = 0; i < 3; i++) {
    chars[i] = getHexCharacter(data & 0x0F);
    data >>= 4;
  }
  str += getHexCharacter(data & 0x0F);
  for (int8_t i = 2; i >= 0; i--) {
    str += chars[i];
  }
  return str;
}

void delay_ms(const uint32_t time) {
  using namespace std::chrono;
  using namespace std::this_thread;

  sleep_for(milliseconds(time));
}
