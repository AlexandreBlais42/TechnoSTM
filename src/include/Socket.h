#ifndef SOCKET_H
#define SOCKET_H

#include <iostream>
#include <cstdint>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <vector>

#include "Error.h"

class SocketServer{
public:
  uint16_t port;
  int fd; // fd veut dire File Descriptor et est un entier qui représente une connection
  struct sockaddr_in address;
  int opt;
  SocketServer();

  void begin(const uint16_t port);
};

#endif // SOCKET_H
