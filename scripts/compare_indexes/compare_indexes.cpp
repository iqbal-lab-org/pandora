#include <iostream>
#include "index.h"

int main(int argc, char* argv[]) {
	Index index1;
	Index index2;
	index1.load(argv[1]);
	index2.load(argv[2]);
	if (index1 == index2) {
		std::cout << "Indexes are EQUAL" << std::endl;
		return 0;
	}
	else {
		std::cout << "Indexes are DIFFERENT" << std::endl;
		return 1;
	}
}
