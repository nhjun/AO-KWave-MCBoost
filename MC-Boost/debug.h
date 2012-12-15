/*
 *  debug.h
 *  MonteCarloThread
 */

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

//#define DEBUG 1

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)
