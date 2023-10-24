#ifndef SOCKET_H
#define SOCKET_H

#include <iostream>
#include <cstdint>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>

class SocketServer{
public:
  uint16_t port;
  int fd, sock;
  struct sockaddr_in address;
  int opt;
  SocketServer();

  void begin(const uint16_t port);
};

#endif // SOCKET_H
