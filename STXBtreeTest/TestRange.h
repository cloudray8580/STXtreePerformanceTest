#ifndef TestRange_HPP
#define TestRange_HPP

// remember to turn off precompiled header in C/C++ precompiled header
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <stx/btree.h>

#include "HilbertDecomposition.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgid = boost::geometry::index::detail;

using namespace std;

typedef bg::model::d2::point_xy<double> MyPoint; // 2d
typedef bg::model::box<MyPoint> Box;

void TestRangeWithRtree(string datafile = "C:/Users/Cloud/Desktop/LearnIndex/data/HilbertSortedPOIs2.csv");
void inline CalculateRangeR(bgi::rtree<MyPoint, bgi::linear<16>>& rtree, string testfile = "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection1km.csv");
void TestRangeWithStx(unsigned max_bits = 14, string datafile = "C:/Users/Cloud/Desktop/LearnIndex/data/HilbertSortedPOIs2.csv");
void inline CalculateRange(stx::btree<long, vector<double>>& btree, unsigned max_bits, string testfile = "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection1km.csv");
void inline RetrieveItems(vector<interval>& intervals, stx::btree<long, vector<double>>& btree, clock_t& time_spent);
vector<vector<double>> ReadSearchAreaFromFile(string inputfile);
vector<vector<double>> LoadDataset(string inputfile);

void TestRangeWithRtree(string datafile) {
	// Load Dataset
	vector<vector<double>> original_dataset = LoadDataset(datafile);

	// Store into boost Rtree
	bgi::rtree<MyPoint, bgi::linear<16>> rtree;
	for (int i = 0; i < original_dataset.size(); i++) {
		MyPoint point(original_dataset[i][0], original_dataset[i][1]);
		rtree.insert(point);
	}

	// query dataset
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection100m.csv");
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection200m.csv");
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection500m.csv");

	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection1km.csv");
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection2km.csv");
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection5km.csv");
	CalculateRangeR(rtree, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection10km.csv");
}

void inline CalculateRangeR(bgi::rtree<MyPoint, bgi::linear<16>>& rtree, string testfile) {
	// Load Testset
	vector<vector<double>> testset = ReadSearchAreaFromFile(testfile);
	clock_t total_time = 0, each_time = 0, start, stop;
	std::vector<MyPoint> returned_values;
	for (int i = 0; i < testset.size(); i++) {
		start = clock();
		returned_values.clear();
		// query range
		Box query_region(MyPoint(testset[i][0], testset[i][2]), MyPoint(testset[i][1], testset[i][3]));
		rtree.query(bgi::intersects(query_region), std::back_inserter(returned_values));
		stop = clock();
		each_time = stop - start;
		total_time += each_time;
	}
	cout << "total time: " << total_time << endl;
}

void TestRangeWithStx(unsigned max_bits, string datafile) {
	// Load Dataset
	vector<vector<double>> original_dataset = LoadDataset(datafile);

	vector<pair<long, vector<double>>> dataset;
	for (int i = 0; i < original_dataset.size(); i++) {
		dataset.push_back(pair<long, vector<double>>(original_dataset[i][4], original_dataset[i]));
	}

	stx::btree<long, vector<double>> btree;
	btree.bulk_load(dataset.begin(), dataset.end());

	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection100m.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection200m.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection500m.csv");

	CalculateRange(btree, 14); // 1km

	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection2km.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection5km.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection10km.csv");
	/*CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection20km.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection50km.csv");
	CalculateRange(btree, 14, "C:/Users/Cloud/Desktop/LearnIndex/data/RangeQueryCollection100km.csv");*/
}



// the max_bits is used to defined the Hilbert sapce
void inline CalculateRange(stx::btree<long, vector<double>>& btree, unsigned max_bits, string testfile) {

	// Load Testset
	vector<vector<double>> testset = ReadSearchAreaFromFile(testfile);

	// count range search time
	int FULL_BITS = 22;
	int differ_bits = FULL_BITS - max_bits;
	clock_t total_time = 0, each_time = 0;
	for (int i = 0; i < testset.size(); i++) {
		// map the original latitude and longitude to integer
		int lower_x = (int)((2 * testset[i][0] + 180) * 10000) >> differ_bits;
		int upper_x = (int)((2 * testset[i][1] + 180) * 10000) >> differ_bits;
		int lower_y = (int)((testset[i][2] + 180) * 10000) >> differ_bits;
		int upper_y = (int)((testset[i][3] + 180) * 10000) >> differ_bits;
		
		// decompost the range into Hilbert interval
		// unsigned max_bits, unsigned lower_x, unsigned upper_x, unsigned lower_y, unsigned upper_y
		vector<interval> intervals = GetIntervals(max_bits, lower_x, upper_x, lower_y, upper_y);
		RetrieveItems(intervals, btree, each_time);
		total_time += each_time;
	}
	cout << "Total Time(ms): " << total_time << endl;
}

//void inline CalculateRange(double lx, double ux, double ly, double uy, unsigned max_bits) {
//
//	int FULL_BITS = 22;
//	int differ_bits = FULL_BITS - max_bits;
//
//	// map the original latitude and longitude to integer
//	int lower_x = (int)((2 * lx + 180) * 10000) >> differ_bits;
//	int upper_x = (int)((2 * ux + 180) * 10000) >> differ_bits;
//	int lower_y = (int)((ly + 180) * 10000) >> differ_bits;
//	int upper_y = (int)((uy + 180) * 10000) >> differ_bits;
//
//	// decompost the range into Hilbert interval
//	// unsigned max_bits, unsigned lower_x, unsigned upper_x, unsigned lower_y, unsigned upper_y
//	vector<interval> intervals = GetIntervals(max_bits, lower_x, upper_x, lower_y, upper_y);
//}

void inline RetrieveItems(vector<interval>& intervals, stx::btree<long, vector<double>>& btree, clock_t& time_spent) {
	clock_t start, stop;
	start = clock();
	vector<vector<double>> results;
	stx::btree<long, vector<double>>::iterator iter;
	int count = 0;
	for (int i = 0; i < intervals.size(); i++) {
		iter = btree.lower_bound(intervals[i].lower_hilbert_value);
		if (iter != btree.end()) {
			while (iter->first <= intervals[i].upper_hilbert_value) {
				results.push_back(iter->second); // consider should I comment it
				count++;
				iter++;
			}
		}
	}
	stop = clock();
	time_spent = stop - start;
	cout << "count: " << count << "  ";
}

// file type: csv (with header line)
// data type: double (latitude and longitude)
// file format: lower_x,upper_x,lower_y,upper_y
vector<vector<double>> ReadSearchAreaFromFile(string inputfile) {

	ifstream infile;
	infile.open(inputfile);

	vector<vector<double>> ranges;
	vector<double> range;
	string temp;

	getline(infile, temp, '\n'); // ignore the header
	while (getline(infile, temp, ',')) {
		range.clear();
		// lower_x [0]
		range.push_back(stod(temp));
		// upper_x [1]
		getline(infile, temp, ',');
		range.push_back(stod(temp));
		// lower_y [2]
		getline(infile, temp, ',');
		range.push_back(stod(temp));
		// upper_y [3]
		getline(infile, temp, '\n'); //last one
		range.push_back(stod(temp));
		// the entire range
		ranges.push_back(range);
	}

	return ranges;
}

// file type: csv (without header line)
// file format: x[double],y[double],Hilbert_x[long],Hilbert_y[long],Hilbert[long],position[int]
vector<vector<double>> LoadDataset(string inputfile) {
	
	ifstream infile;
	infile.open(inputfile); // typical: "C:/Users/Cloud/Desktop/LearnIndex/data/HilbertSortedPOIs2.csv"
	
	vector<vector<double>> dataset;
	vector<double> data;
	string temp;

	int count = 0;

	while (getline(infile, temp, ',')) {
		data.clear();
		// original_x [0]
		data.push_back(stod(temp));
		// original_y [1]
		getline(infile, temp, ',');
		data.push_back(stod(temp));
		// Hilbert_x [2]
		getline(infile, temp, ',');
		data.push_back(stod(temp));
		// Hilbert_y [3]
		getline(infile, temp, ',');
		data.push_back(stod(temp));
		// Hilbert_value [4]
		getline(infile, temp, ',');
		data.push_back(stod(temp));
		// position [5]
		getline(infile, temp, '\n'); //last one
		data.push_back(stod(temp));

		dataset.push_back(data);
		count++;
		if (count % 10000 == 0) {
			cout << "loading: " << count << endl;
		}
	}

	return dataset;
}

#endif // TestRange_HPP
