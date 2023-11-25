#ifndef MAIN_H 
	#define MAIN_H
	using namespace std;

	#include <iostream> 
	#include <math.h> 
	#include <cstdlib>
	#include <algorithm>
	#include "gurobi_c++.h"
	#include <ilcp/cp.h>
	#include "time.h"
	#include "Allocation.h" 

	float EPSILON = 0.001;
	
	int step0(Allocation& allo, vector<Item>& itemsCSP, vector<int>& demandCSP, vector<vector<int> >& indexCSP);
	int step1(Allocation& allo, vector<vector<vector<int> > > & cuts1, const vector<Item>& itemsCSP, const vector<int>& demandCSP, const vector<vector<int> >& indexCSP);
	int step2(Allocation& allo, const vector<vector<int> > & cut1, vector<vector<vector<int> > > & cuts2, const vector<Item>& itemsCSP, const vector<vector<int> >& indexCSP);
	int step3(Allocation& allo, const vector<vector<int> > & cut2, const vector<Item>& itemsCSP, const vector<vector<int> >& indexCSP);
		
	int preprocessingWDP(vector<Item>& items, int& W);
	int preprocessingLDP(vector<Item>& items, int& L);
	int preprocessingHDP(vector<Item>& items, int& H);	
	int preprocessingWLP(vector<Item>& items, int& W);
	int preprocessingHLP(vector<Item>& items, int& H, const int& maxStack);	
#endif 
