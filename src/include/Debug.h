#ifdef DEBUG_MODE
#define DEBUG std::cout << "[DEBUG] : "
#else
#define DEBUG                                                                  \
  if (false)                                                                   \
  std::cout
#endif // DEBUG_MODE
