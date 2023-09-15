#ifndef SERIAL_H
#define SERIAL_H 1

#include <string>
#include <iostream>

class Serial{
public:
  Serial();

  void begin(const uint32_t baudrate);
  void write(const std::string s);
  
}

#endif //SERIAL_H 
