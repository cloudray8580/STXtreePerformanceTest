// STXBtreeTest.cpp
#ifndef STXBtreeTest_CPP
#define STXBtreeTest_CPP

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <stx/btree.h> // if use this, need to turn the x86 -> x64!!!
//#include "D:\stxbtree\stx-btree-0.9\include\stx\btree.h" // do not need to change to x64

#include "TestRange.h"
#include "TestPoint.h"

using namespace std;

int main()
{
	//TestPointLookup();

	//TestRangeWithRtree();
	TestRangeWithStx();

	system("pause");
}

#endif