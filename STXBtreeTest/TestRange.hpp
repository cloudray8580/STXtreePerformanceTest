
#ifndef TestRange_HPP
#define TestRange_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <stx/btree.h>

#include "HilbertDecomposition.h"

using namespace std;

void TestRangeWithStx() {
	bitmask_t coord[2] = { 2244891,2941279 };
	cout << hilbert_c2i(2, 22, coord) << endl;
}

#endif // TestRange_HPP
