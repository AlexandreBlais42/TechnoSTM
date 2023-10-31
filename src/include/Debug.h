#ifndef DEBUG_H
#define DEBUG_H


#ifdef DEBUG_MODE
#define DEBUG std::cout << "[DEBUG] : "
#else
#define DEBUG                                                                  \
  if (false)                                                                   \
  std::cout
#endif // DEBUG_MODE

#endif // DEBUG_H
