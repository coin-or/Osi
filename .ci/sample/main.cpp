#include <iostream>
#include <coin/CoinUtilsConfig.h>
#include <coin/OsiConfig.h>

int main(int argc, char** argv) {
  std::cout << "coinutils version: " << COINUTILS_VERSION << std::endl;
  std::cout << "osi version: " << OSI_VERSION << std::endl;
	return 0;
}

