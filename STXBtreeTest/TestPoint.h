#pragma once
#ifndef Trash_H
#define Trash_H

#include <iostream>
#include <fstream>

union DatasetItem {
	double dvalue;
	unsigned long ivalue;
};

int TEST_CYCLE = 10;
int TEST_NUM = 10000;

void TestPointLookup() {

	ifstream infile;
	infile.open("C:/Users/Cloud/Desktop/LearnIndex/data/HilbertSortedPOIs2.csv");
	int countline = 0;
	string latitude_s;
	string longitude_s;
	string latitude_h_s;
	string longitude_h_s;
	string hilbert_s;
	string order_s;

	/*double latitude;
	double longitude;
	unsigned long latitude_h;
	unsigned long longitude_h;
	unsigned long hilbert;
	unsigned long order;*/

	DatasetItem latitude;
	DatasetItem longitude;
	DatasetItem latitude_h;
	DatasetItem longitude_h;
	DatasetItem hilbert;
	DatasetItem order;

	vector<pair<long, vector<DatasetItem>>> dataset;
	vector<DatasetItem> onerow;

	cout << "data loading..." << endl;
	bool firstFlag = true;
	while (getline(infile, latitude_s, ',')) {
		onerow.clear();

		latitude.dvalue = stod(latitude_s);
		onerow.push_back(latitude);

		getline(infile, longitude_s, ',');
		longitude.dvalue = stod(longitude_s);
		onerow.push_back(longitude);

		getline(infile, latitude_h_s, ',');
		latitude_h.ivalue = stol(latitude_h_s);
		onerow.push_back(latitude_h);

		getline(infile, longitude_h_s, ',');
		longitude_h.ivalue = stol(longitude_h_s);
		onerow.push_back(longitude_h);

		getline(infile, hilbert_s, ',');
		hilbert.ivalue = stol(hilbert_s);
		onerow.push_back(hilbert);

		getline(infile, order_s, '\n'); //last one
		order.ivalue = stol(order_s);
		onerow.push_back(order);


		dataset.push_back(pair<long, vector<DatasetItem>>(hilbert.ivalue, onerow));
		cout << std::fixed;
		if (countline % 10000 == 0) {
			if (!firstFlag) {
				cout << "\r";
			}
			cout << setprecision(2) << (double)countline * 100 / 1048575 << "%";
			firstFlag = false;
		}

		/*if (countline >= 1048573) {
			cout << "number of lines: " << countline << endl;
			cout << "latitude: " << latitude_s << endl;
			cout << "longitude: " << longitude_s << endl;
			cout << "latitude_h: " << latitude_h_s << endl;
			cout << "longitude_h: " << longitude_h_s << endl;
			cout << "hilbert: " << hilbert_s << endl;
			cout << "order: " << order_s << endl;
		}*/

		/*if (countline < 10) {
			cout << "number of lines: " << countline << endl;
			cout << "latitude: " << latitude << endl;
			cout << "longitude: " << longitude << endl;
			cout << "latitude_h: " << latitude_h << endl;
			cout << "longitude_h: " << longitude_h << endl;
			cout << "hilbert: " << hilbert << endl;
			cout << "order: " << order << endl;
		}
		else {
			break;
		}*/

		countline++;
	}
	cout << endl << "data loading finished..." << endl;
	stx::btree<long, vector<DatasetItem>> btree;
	btree.bulk_load(dataset.begin(), dataset.end());
	//cout << "data stored in stx btree finished..." << endl;

	int data_size = dataset.size();
	map<long, vector<DatasetItem>> stlmap;
	for (int i = 0; i < data_size; i++) {
		stlmap.insert(dataset[i]);
	}

	vector<long> testset;
	srand((unsigned)time(NULL));
	clock_t start1, stop1, total1 = 0; // CLOCKS_PER_SEC = 1000
	clock_t start2, stop2, total2 = 0; // CLOCKS_PER_SEC = 1000
	for (int test_cycle = 0; test_cycle < TEST_CYCLE; test_cycle++) {
		testset.clear();
		long index = 0;
		for (int i = 0; i < TEST_NUM; i++) {
			index = rand() % data_size;
			testset.push_back(dataset[index].first);
		}

		// Query And Time Counting
		stx::btree<long, vector<DatasetItem>>::iterator iter;

		start1 = clock();
		auto t0 = chrono::steady_clock::now();

		//iter = btree.find(testset[0]); // average 300 - 900 ns
		for (int i = 0; i < TEST_NUM; i++) { // average 50k ns
			iter = btree.find(testset[i]);
		}

		auto t1 = chrono::steady_clock::now();
		stop1 = clock();
		
		//total1 += stop1 - start1;
		cout << "# Query Case = " << test_cycle << " Total query time for stx btree: " << stop1 - start1 << endl;
		cout << "# Query Case = " << test_cycle << ": Total Time in chrono For STXtree: " << chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() << " ns" << endl;

		map<long, vector<DatasetItem>>::iterator iter2;

		start2 = clock();
		t0 = chrono::steady_clock::now();

		//iter2 = stlmap.find(testset[0]);
		for (int i = 0; i < TEST_NUM; i++) {
			iter2 = stlmap.find(testset[i]);
		}

		t1 = chrono::steady_clock::now();
		stop2 = clock();
		
		//total2 += stop2 - start2;
		cout << "# Query Case = " << test_cycle << "Total query time for STL map: " << stop2 - start2 << endl;
		cout << "# Query Case = " << test_cycle << ": Total Time in chrono For STL map: " << chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count() << " ns" << endl;
	}
	//cout << "average of above stx test cases: " << total1 / TEST_CYCLE << endl;
	//cout << "average of above STL test cases: " << total2 / TEST_CYCLE << endl;

}
#endif
