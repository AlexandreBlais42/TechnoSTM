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
  int fd;
  struct sockaddr_in address;
  int opt;
  std::vector<int> connections();
  SocketServer();

  void begin(const uint16_t port);
  void acceptConnections();
};

#endif // SOCKET_H
