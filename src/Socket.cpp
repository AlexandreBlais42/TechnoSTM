#include "Socket.h"
#include <sys/socket.h>

SocketServer::SocketServer(){}

void SocketServer::begin(const uint16_t port){
  this->port = port;
  if ((fd = socket(AF_INET, SOCK_STREAM, 0)) < 0){
    COUTERROR;
    COUTERRNO;
    exit(1);
  }

  if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt))){
    std::cout << "Erreur lors de l'appel de setsockopt()" << std::endl;
  }
  
  address.sin_family = AF_INET;
  address.sin_addr.s_addr = INADDR_ANY;
  address.sin_port = htons(port);

  if (bind(fd, (struct sockaddr *) &address, sizeof(address)) < 0){
    COUTERROR;
    COUTERRNO;
  }

  if(listen(fd, 1) < 0){ // N'autorise qu'une seule connection en backlog
    COUTERROR;
    COUTERRNO;
  }
}
