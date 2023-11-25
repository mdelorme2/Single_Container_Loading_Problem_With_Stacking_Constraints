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

	class mycallbackD2: public GRBCallback
	{
		public:
			int ncuts;
			Allocation allo;
			vector<GRBVar> isItemPackedV;
			vector<vector<vector<GRBVar> > > isItemCSPPackedInCoN;
			vector<vector<vector<GRBVar> > > isItemCSPPackedInCoR;
			vector<vector<int> > indexCSP;
			int maxV;
			
			mycallbackD2(const Allocation& xallo, const vector<vector<int> >& xindexCSP, const vector<GRBVar>& xisItemPackedV, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoN, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoR);
			int getNcuts();
			int getMaxV();
			vector<vector<int> > getPps();
			int bestSol;

		protected:
			void callback();
	};

	
	int model(Allocation& allo);
	int step2(Allocation& allo, const vector<vector<int> >& indexCSP, const vector<vector<int> >& cut);
	
	int preprocessingWDP(vector<Item>& items, int& W);
	int preprocessingLDP(vector<Item>& items, int& L);
	int preprocessingHDP(vector<Item>& items, int& H);	
	int preprocessingWLP(vector<Item>& items, int& W);
	int preprocessingHLP(vector<Item>& items, int& H, const int& maxStack);	
#endif 
