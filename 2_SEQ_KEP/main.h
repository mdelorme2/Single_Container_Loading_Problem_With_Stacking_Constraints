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
	#include <ilcp/cp.h>
	
	float EPSILON = 0.001;

	class mycallbackKEP: public GRBCallback
	{
		public:
			int ncuts;
			Allocation allo;
			vector<vector<GRBVar> > isArcUsed;
			vector<GRBVar> isItemPackedV;
			
			mycallbackKEP(const Allocation& xallo, const vector<GRBVar>& xisItemPackedV, const vector<vector<GRBVar> >& xisArcUsed);
			int getNcuts();

		protected:
			void callback();
	};

	class mycallbackYC: public GRBCallback
	{
		public:
			int ncuts;
			Allocation allo;
			vector<vector<vector<GRBVar> > > isItemCSPPackedInCoN, isItemCSPPackedInCoR;
			vector<Item> itemsCSP;
			vector<vector<vector<int> > >ordinates;
			vector<int> NPW;
			mycallbackYC(const Allocation& xallo, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoN, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoR, const vector<Item>& xitemsCSP, const vector<int>& xNPW);
			int getNcuts();
			int bestSol;
			vector<vector<vector<int> > >getOrdinates();
		protected:
			void callback();
	};
	
	int step1(Allocation& allo);
	int step2(Allocation& allo);
		
	int preprocessingWDP(vector<Item>& items, int& W);
	int preprocessingLDP(vector<Item>& items, int& L);
	int preprocessingHDP(vector<Item>& items, int& H);	
	
	int preprocessingWLP(vector<Item>& items, int& W);
	int preprocessingHLP(vector<Item>& items, int& H, const int& maxStack);	
#endif 
