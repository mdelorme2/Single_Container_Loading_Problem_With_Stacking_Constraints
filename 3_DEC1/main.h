#ifndef MAIN_H 
	#define MAIN_H
	using namespace std;

	#include <iostream> 
	#include <math.h> 
	#include <cstdlib>
	#include <algorithm>
	#include "gurobi_c++.h"
	#include "time.h"
	#include "Allocation.h" 

	float EPSILON = 0.001;
	
	class mycallbackD1: public GRBCallback
	{
		public:
			int ncuts;
			Allocation allo;
			vector<GRBVar> isItemPackedV;
			vector<vector<GRBVar> > isItemCSPPackedV;
			int maxV;
			
			mycallbackD1(const Allocation& xallo, const vector<GRBVar>& xisItemPackedV, const vector<vector<GRBVar> >& xisItemCSPPackedV);
			int getNcuts();
			int getMaxV();
			vector<vector<int> > getPps();
			int bestSol;

		protected:
			void callback();
	};


	int model(Allocation& allo);
	int preprocessingWDP(vector<Item>& items, int& W);
	int preprocessingLDP(vector<Item>& items, int& L);
	int preprocessingHDP(vector<Item>& items, int& H);	
	int preprocessingWLP(vector<Item>& items, int& W);
	int preprocessingHLP(vector<Item>& items, int& H, const int& maxStack);	
	int step2(Allocation& allo);
#endif 
