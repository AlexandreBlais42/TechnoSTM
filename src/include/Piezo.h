#ifndef PIEZO_H
#define PIEZO_H 1

#include <iostream>
#include <cstdint>

class Piezo{
public:
  Piezo();
  void setValue(const uint16_t value);
};

#endif //PIEZO_H
