#include "Utils.h"

char map4BitToHexChar(const uint8_t val){
  assert(val < (1 << 4));
  if (val <= 9){
    return '0' + val;
  }
  return 'A' + val - 10;
}

std::string getHexString(uint16_t data){
  std::string str = "";
  char chars[3];
  for (uint8_t i = 0 ; i < 3 ; i++){
    chars[i] = map4BitToHexChar(data & 0x0F);
    data >>= 4;
  }
  str += map4BitToHexChar(data & 0x0F);
  for (int8_t i = 3 ; i >= 0 ; i--){
    str += chars[i];
  }
  return str;
}

void delay_ms(const uint32_t time){
  using namespace std::chrono;
  using namespace std::this_thread;

  sleep_for(milliseconds(time));
}
