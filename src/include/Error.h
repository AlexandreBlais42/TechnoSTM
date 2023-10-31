
#ifndef ERROR_H
#define ERROR_H

#include <cstring>

#define COUTERROR std::cout << "\033[31m[ERROR] " << __PRETTY_FUNCTION__ << " près de la ligne " << __LINE__ << "\033[0m:" << std::endl
#define COUTERRNO std::cout << "\033[31merrno: " << strerror(errno) << "\033[0m" << std::endl

#endif // ERROR_H
