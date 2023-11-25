#include "main.h"

/*	*************************************************************************************
	*************************************  MAIN *****************************************
	************************************************************************************* */

double initTimeModelCPU;

int main(int argc, char **argv){
	
	initTimeModelCPU = getCPUTime();
	
	// local variables
	Allocation allo ;
	string filein = argv[2];
	string path = argv[1];
	string pathAndFileout = argv[3];
	string picture = argv[4];

	// functions
	allo.load(path,filein);
	allo.printProb();
	preprocessingWLP(allo.items, allo.W);
	preprocessingLDP(allo.items, allo.L);
	preprocessingHLP(allo.items, allo.H, allo.maxStack);
	step1(allo);
	step2(allo);
	allo.rescaleW();
	allo.printInfo(pathAndFileout);
	allo.printFile(picture);
}

int step1(Allocation& allo){
	
	// Transform into CSP
	vector<Item> itemsCSP;
	vector<int> demandCSP;
	vector<vector<int> > indexCSP;
	for(int i =0; i<allo.items.size();i++){
		bool exist = false;
		for(int j = 0; j<itemsCSP.size();j++){
			if(allo.items[i].len == itemsCSP[j].len && allo.items[i].wid == itemsCSP[j].wid && allo.items[i].lenR == itemsCSP[j].lenR && allo.items[i].widR == itemsCSP[j].widR){
				cout << "item " << i  << " is the same as item " << indexCSP[j][0] << endl;
				exist = true;
				demandCSP[j] += 1;
				indexCSP[j].push_back(i);
				break;
			}
		}
		if(!exist){
			itemsCSP.push_back(allo.items[i]);
			demandCSP.push_back(1);
			indexCSP.push_back({i});
		}
	}	

	for (int i = 0; i < indexCSP.size(); i++){	
		for (int j = 0; j < indexCSP[i].size(); j++){	
			cout << indexCSP[i][j] << " ";
		}
		cout << endl;
	}
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 16);
	env.start();
		
	// Step 1 -- make columns
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;

		// Related to the item type 
		vector<GRBVar> isItemCSPPackedV(itemsCSP.size());		
		vector<GRBLinExpr> isItemPacked(itemsCSP.size(),0);
		
		// Related to the item 
		vector<GRBVar> isItemPackedV(allo.nbItems);
		vector<GRBVar> isItemStackedV(allo.nbItems);
		vector<vector<GRBVar> > isArcUsed (allo.cArcs.size());
		vector<vector<GRBLinExpr> > cIn (allo.nbItems);
		vector<vector<GRBLinExpr> > cOut (allo.nbItems);
		vector<GRBLinExpr> isItemStacked(allo.nbItems,0);

		// Initialization
		for (int i = 0; i < itemsCSP.size(); i++)
			isItemCSPPackedV[i] = model.addVar(0, demandCSP[i], 0, GRB_INTEGER);
		
		for (int i = 0; i < allo.items.size(); i++){	
			isItemStackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			cIn[i].resize(allo.maxStack,0);
			cOut[i].resize(allo.maxStack,0);
		}

		for (int i = 0; i < allo.cArcs.size(); i++){	
			isArcUsed[i].resize(allo.maxStack-1);
			for (int j = 0; j < allo.maxStack-1; j++)
				isArcUsed[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
		}
		
		model.update();

		// Compute values packing
		for (int i = 0; i < itemsCSP.size(); i++){			
			for(int j = 0; j < indexCSP[i].size();j++)
				isItemPacked[i] += isItemPackedV[indexCSP[i][j]];
		}
		
		//  Compute values stacking 
		for (int i = 0; i < allo.cArcs.size(); i++){
			for (int j = 0; j < allo.maxStack-1; j++){
				cOut[allo.cArcs[i][0]][j] += isArcUsed[i][j];
				cIn[allo.cArcs[i][1]][j+1] += isArcUsed[i][j];
				isItemStacked[allo.cArcs[i][1]] += isArcUsed[i][j];
			}
		}

		// Global constraints -- linking var between stacking and packing
		for (int i = 0; i < itemsCSP.size(); i++)
			model.addConstr(isItemPacked[i] == isItemCSPPackedV[i]);
		
		//  Global constraints -- stacking
		for (int i = 0; i < allo.items.size(); i++){
			model.addConstr(cOut[i][0] <= isItemPackedV[i]);
			for (int j = 1; j < allo.maxStack-1; j++)
				model.addConstr(cOut[i][j] <= cIn[i][j]);
			model.addConstr(isItemStacked[i] == isItemStackedV[i]);
		}

		// Global constraints -- one item exactly
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(isItemPackedV[i] + isItemStackedV[i] == 1);

		// Objective function -- minimize column area
		for (int i = 0; i < itemsCSP.size(); i++)
			objFun += isItemCSPPackedV[i] * (itemsCSP[i].lenO * itemsCSP[i].widO);
		
		model.setObjective(objFun, GRB_MINIMIZE);
			
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
		
		// Callback
		mycallbackKEP cb = mycallbackKEP(allo, isItemPackedV, isArcUsed);
		model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
		model.setCallback(&cb);                     		// link the callback to the model	
		model.optimize();
		allo.infos.nbCuts1 = cb.getNcuts();

		// Identify the columns
		vector<int> succ (allo.nbItems,-1);
		for (int i = 0; i < allo.cArcs.size(); i++){
			for (int j = 0; j < allo.maxStack-1; j++){
				if(ceil(isArcUsed[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1)
					succ[allo.cArcs[i][0]] = allo.cArcs[i][1];
			}
		}
				
		// Filling Info and solution
		allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
		allo.columns.resize(allo.items.size());
		for (int i = 0; i < allo.nbItems; i++){
			if(ceil(isItemPackedV[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				allo.columns[i].push_back(i); 
				int curr = i; 
				while(succ[curr] != -1){ 
					allo.columns[i].push_back(succ[curr]); 
					curr = succ[curr];
				}
			}
		}		
	}
	
	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}	
		
	// End
	return 0;
}
		
int step2(Allocation& allo){
	// Create reduced instance
	vector<Item> itemsBottom;  
	vector<int> values;

	int value = 0;
	for (int i = 0; i < allo.columns.size(); i++){
		if(allo.columns[i].size() > 0){
			cout << "Column " << i << ": [";
			itemsBottom.push_back(allo.items[allo.columns[i][0]]);
			int value =  0;
			for(int j = 0; j < allo.columns[i].size();j++){
				cout << allo.columns[i][j] << " ";
				value += allo.items[allo.columns[i][j]].pro;
			}
			values.push_back(value);
			cout << "]" << "(" << value << ") " << endl;
		}
	}
	cout << endl;

	preprocessingWDP(itemsBottom, allo.W);
	preprocessingLDP(itemsBottom, allo.L);	
	
	// Transform into CSP
	vector<Item> itemsCSP;
	vector<int> demandCSP;
	vector<vector<int> > indexCSP;
	for(int i =0; i<itemsBottom.size();i++){
		bool exist = false;
		for(int j = 0; j<itemsCSP.size();j++){
			if(itemsBottom[i].len == itemsCSP[j].len && itemsBottom[i].wid == itemsCSP[j].wid && itemsBottom[i].lenR == itemsCSP[j].lenR && itemsBottom[i].widR == itemsCSP[j].widR){
				cout << "itemsBottom " << i  << " is the same as item " << indexCSP[j][0] << endl;
				exist = true;
				demandCSP[j] += 1;
				indexCSP[j].push_back(i);
				break;
			}
		}
		if(!exist){
			itemsCSP.push_back(itemsBottom[i]);
			demandCSP.push_back(1);
			indexCSP.push_back({i});
		}
	}	

	for (int i = 0; i < indexCSP.size(); i++){	
		for (int j = 0; j < indexCSP[i].size(); j++){	
			cout << indexCSP[i][j] << " ";
		}
		cout << endl;
	}
		
	// Local variables
	vector<int> NPW (allo.W+1,0);
	int minW = allo.W;
	
	// Compute normal patterns
	NPW[0] = 1;
	for(int i =0; i<itemsCSP.size();i++){
		for(int j = 0; j < demandCSP[i]; j++){
			minW = min(minW,itemsCSP[i].wid); minW = min(minW,itemsCSP[i].widR); 
			for(int k = allo.W - min(itemsCSP[i].wid,itemsCSP[i].widR); k>=0; k--){
				if(NPW[k] == 1 && k + itemsCSP[i].wid <= allo.W) NPW[k + itemsCSP[i].wid] = 1;
				if(NPW[k] == 1 && k + itemsCSP[i].widR <= allo.W) NPW[k + itemsCSP[i].widR] = 1;
			}
		}
	}

	// Remove ending NPs
	for(int i = allo.W - minW + 1; i < allo.W; i++) NPW[i] = 0;
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 16);
	env.start();

	// Step 1 -- make columns
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		vector<GRBLinExpr> columnFilled (allo.W,0);
		GRBLinExpr objFun = 0;

		// Related to the item
		vector<GRBVar> isColumnSeclected(itemsBottom.size());
		vector<GRBLinExpr> isColumnSeclectedT(itemsCSP.size(),0);
				
		// Related to the item type 
		vector<GRBVar> isItemCSPPackedV(itemsCSP.size());			
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoN (itemsCSP.size());
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoR (itemsCSP.size());
		vector<GRBLinExpr> isItemCSPPacked(itemsCSP.size(),0);		

		// Initialization
		for (int i = 0; i < itemsCSP.size(); i++){
			isItemCSPPackedV[i] = model.addVar(0, demandCSP[i], 0, GRB_INTEGER);
			isItemCSPPackedInCoN[i].resize(ceil(log2(demandCSP[i]))+1);
			isItemCSPPackedInCoR[i].resize(ceil(log2(demandCSP[i]))+1);
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				isItemCSPPackedInCoN[i][j].resize(allo.W);
				isItemCSPPackedInCoR[i][j].resize(allo.W);
				for (int k = 0; k < allo.W; k++){
					if(NPW[k] && k + itemsCSP[i].wid <= allo.W){
						if(k + itemsCSP[i].wid != allo.W && k + itemsCSP[i].wid + minW > allo.W && NPW[allo.W - itemsCSP[i].wid])
							isItemCSPPackedInCoN[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
						else
							isItemCSPPackedInCoN[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);
					}
					if(NPW[k] && k + itemsCSP[i].widR <= allo.W){
						if(itemsCSP[i].lenR == itemsCSP[i].len && itemsCSP[i].widR == itemsCSP[i].wid)
							isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
						else{
							if(k + itemsCSP[i].widR != allo.W && k + itemsCSP[i].widR + minW > allo.W && NPW[allo.W - itemsCSP[i].widR])
								isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
							else
								isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);
						}								
					}				
				}
			}
		}
		for (int i = 0; i < itemsBottom.size(); i++) isColumnSeclected[i]  =  model.addVar(0, 1, 0, GRB_BINARY);
		
		model.update();

		// Compute values packing
		for (int i = 0; i < itemsCSP.size(); i++){			
			if(demandCSP[i] == 0) continue;		
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				for (int k = 0; k < allo.W; k++){
					if(NPW[k] && k + itemsCSP[i].wid <= allo.W){
						isItemCSPPacked[i] += pow(2,j) * isItemCSPPackedInCoN[i][j][k];
						for(int k2 = 0; k2 < itemsCSP[i].wid ; k2++){
							if(NPW[k+k2])
								columnFilled[k+k2] += pow(2,j) * itemsCSP[i].len * isItemCSPPackedInCoN[i][j][k];
						}
					}
					if(NPW[k] && k + itemsCSP[i].widR <= allo.W){
						isItemCSPPacked[i] += pow(2,j) * isItemCSPPackedInCoR[i][j][k];
						for(int k2 = 0; k2 < itemsCSP[i].widR ; k2++){
							if(NPW[k+k2])
								columnFilled[k+k2] += pow(2,j) * itemsCSP[i].lenR * isItemCSPPackedInCoR[i][j][k];
						}
					}
				}
			}
			for(int j = 0; j < indexCSP[i].size();j++)
				isColumnSeclectedT[i] += isColumnSeclected[indexCSP[i][j]];
			model.addConstr(isItemCSPPacked[i] == isItemCSPPackedV[i]);
		}

		// Global constraints -- match column type and column
		for (int i = 0; i < itemsCSP.size(); i++)
			model.addConstr(isItemCSPPackedV[i] == isColumnSeclectedT[i]);
		
		// Global constraints -- length
		for (int i = 0; i < allo.W; i++){
			if(NPW[i])
				model.addConstr(columnFilled[i] <= allo.L);
		}

		// Objective function -- maximize value of selected columns
		for (int i = 0; i < itemsBottom.size(); i++) 
			objFun += values[i] * isColumnSeclected[i];
		
		model.setObjective(objFun, GRB_MAXIMIZE);
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600 - (getCPUTime() - initTimeModelCPU));
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);

		// Callback
		mycallbackYC cb = mycallbackYC(allo, isItemCSPPackedInCoN, isItemCSPPackedInCoR, itemsCSP, NPW);
		model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
		model.setCallback(&cb);                     		// link the callback to the model	
		model.optimize();
		allo.infos.nbCuts2 = cb.getNcuts();
		vector<vector<vector<int> > > ordinates = cb.getOrdinates();
		
		// Filling Info
		allo.infos.UB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		allo.infos.opt = false;
		allo.infos.LB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);	
		if(allo.infos.LB == allo.infos.UB) allo.infos.opt = true;
		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);		
		allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
		
		// Filling Solution
		allo.pps.resize(allo.items.size());
		for (int i = 0; i < itemsCSP.size(); i++){	
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				for (int k = 0; k < allo.W; k++){
					if(NPW[k] && k + itemsCSP[i].wid <= allo.W && ceil(isItemCSPPackedInCoN[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						for (int l = 0; l < pow(2,j); l++){
							int index = indexCSP[i].back();
							while(ceil(isColumnSeclected[index].get(GRB_DoubleAttr_X) - EPSILON) == 0){
								indexCSP[i].pop_back();
								index = indexCSP[i].back();
							}
							// get ordinates
							for(int l = 0; l < ordinates[i].size();l++){
								if(ordinates[i][l][1] == k && ordinates[i][l][2] == 0){
									allo.pps[itemsBottom[index].idx].push_back(ordinates[i][l][0]);
									ordinates[i][l][2] = -1;
									break;
								}
							}
							allo.pps[itemsBottom[index].idx].push_back(k); allo.pps[itemsBottom[index].idx].push_back(0); 
							indexCSP[i].pop_back();
						}
					}
				}
			}
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				for (int k = 0; k < allo.W; k++){
					if(NPW[k] && k + itemsCSP[i].widR <= allo.W && ceil(isItemCSPPackedInCoR[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						for (int l = 0; l < pow(2,j); l++){
							int index = indexCSP[i].back();	
							while(ceil(isColumnSeclected[index].get(GRB_DoubleAttr_X) - EPSILON) == 0){
								indexCSP[i].pop_back();
								index = indexCSP[i].back();
							}
							// get ordinates
							for(int l = 0; l < ordinates[i].size();l++){
								if(ordinates[i][l][1] == k && ordinates[i][l][2] == 1){
									allo.pps[itemsBottom[index].idx].push_back(ordinates[i][l][0]);
									ordinates[i][l][2] = -1;
									break;
								}
							}
							allo.pps[itemsBottom[index].idx].push_back(k); allo.pps[itemsBottom[index].idx].push_back(1); 
							indexCSP[i].pop_back();
						}
					}
				}
			}
		}
	}
	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}	
		
	// End
	return 0;
}	

int preprocessingWDP(vector<Item>& items, int& W){
	int nbItems = items.size();
	// Preprocessing 1 -- decrease W
	vector<bool> isR(W+1,false); isR[0] = true;
	for (int i = 0; i < nbItems; i++) {
		for(int k = W - min(items[i].widR,items[i].wid); k>=0; k--){
			if(isR[k] == 1 && k + items[i].widR <= W) isR[k + items[i].widR] = 1;
			if(isR[k] == 1 && k + items[i].wid <= W) isR[k + items[i].wid] = 1;
		}
	}
	cout << "W from " << W << " to ";
	while (!isR[W]) W--;
	cout << W << endl;
	
	// Preprocessing 2 -- increase wid
	for (int i = 0; i < nbItems; i++){
		isR.resize(0); isR.resize(W+1,false); isR[0] = true;
		for (int j = 0; j < nbItems; j++){
			if(j != i){
				for(int k = W - min(items[j].wid,items[j].widR); k>=0; k--){
					if(isR[k] == 1 && k + items[j].wid <= W) isR[k + items[j].wid] = 1;
					if(isR[k] == 1 && k + items[j].widR <= W) isR[k + items[j].widR] = 1;
				}
			}
		}
		cout << "Item " << i << " wid from " << items[i].wid << " to ";
		while(!isR[W - items[i].wid]) items[i].wid++;
		cout << items[i].wid << endl;
		cout << "Item " << i << " widR from " << items[i].widR << " to ";
		while(!isR[W - items[i].widR]) items[i].widR++;
		cout << items[i].widR << endl;
	}
	return 0;
}

int preprocessingLDP(vector<Item>& items, int& L){
	int nbItems = items.size();
	// Preprocessing 1 -- decrease L
	vector<bool> isR(L+1,false); isR[0] = true;
	for (int i = 0; i < nbItems; i++) {
		for(int k = L - min(items[i].len,items[i].lenR); k>=0; k--){
			if(isR[k] == 1 && k + items[i].len <= L) isR[k + items[i].len] = 1;
			if(isR[k] == 1 && k + items[i].lenR <= L) isR[k + items[i].lenR] = 1;
		}
	}
	cout << "L from " << L << " to ";
	while (!isR[L]) L--;
	cout << L << endl;
	
	// Preprocessing 2 -- increase length
	for (int i = 0; i < nbItems; i++){
		isR.resize(0); isR.resize(L+1,false); isR[0] = true;
		for (int j = 0; j < nbItems; j++){
			if(j != i){
				for(int k = L - min(items[j].len,items[j].lenR); k>=0; k--){
					if(isR[k] == 1 && k + items[j].len <= L) isR[k + items[j].len] = 1;
					if(isR[k] == 1 && k + items[j].lenR <= L) isR[k + items[j].lenR] = 1;
				}
			}
		}
		cout << "Item " << i << " len from " << items[i].len << " to ";
		while(!isR[L - items[i].len]) items[i].len++;
		cout << items[i].len << endl;
		cout << "Item " << i << " lenR from " << items[i].lenR << " to ";
		while(!isR[L - items[i].lenR]) items[i].lenR++;
		cout << items[i].lenR << endl;
	}
	return 0;
}

int preprocessingHDP(vector<Item>& items, int& H){
	int nbItems = items.size();
	// Preprocessing 1 -- decrease H
	vector<bool> isR(H+1,false); isR[0] = true;
	for (int i = 0; i < nbItems; i++) {
		for(int k = H - items[i].hei; k>=0; k--)
			if(isR[k] == 1) isR[k + items[i].hei] = 1;
	}
	cout << "H from " << H << " to ";
	while (!isR[H]) H--;
	cout << H << endl;
	
	// Preprocessing 2 -- increase height
	for (int i = 0; i < nbItems; i++){
		isR.resize(0); isR.resize(H+1,false); isR[0] = true;
		for (int j = 0; j < nbItems; j++){
			if(j != i){
				for(int k = H - items[j].hei; k>=0; k--)
					if(isR[k] == 1) isR[k + items[j].hei] = 1;
			}
		}
		cout << "Item " << i << " height from " << items[i].hei << " to ";
		while(!isR[H - items[i].hei]) items[i].hei++;
		cout << items[i].hei << endl;
	}
	return 0;
}

int preprocessingWLP(vector<Item>& items, int& W){
	int nbItems = items.size();
	// Local variables
	vector<vector<int> > patternsF; patternsF.push_back({});
	vector<int> patternsC; patternsC.push_back(0);
	for (int i = 0; i < items.size(); i++){	
		vector<vector<int> > patternsToAdd;
		for(int j = 0; j < patternsF.size();j++){
			if(patternsC[j] + items[i].lenO <= W){
				patternsToAdd.push_back(patternsF[j]); patternsToAdd.back().push_back(i);
				patternsC.push_back(patternsC[j] + items[i].lenO);
			}
			if(patternsC[j] + items[i].widO <= W){
				patternsToAdd.push_back(patternsF[j]); patternsToAdd.back().push_back(i+nbItems);
				patternsC.push_back(patternsC[j] + items[i].widO);
			}
		}
		for(int j = 0; j < patternsToAdd.size(); j++){
			patternsF.push_back(patternsToAdd[j]);
		}
	}
	
	/*cout << "Print patterns" << endl;
	for(int i = 0; i < patternsF.size(); i++){
		cout << i << " " << patternsC[i] << "-";
		for(int j = 0; j < patternsF[i].size();j++){
			cout << patternsF[i][j] << " ";
		}
		cout << endl;
	}*/
	cout << "There are " << patternsF.size() << " feasible patterns " << endl;
	if(W <= 50){
		cout << "W is below 50, not much is to be gained by applying the LP preprocessing" << endl;
		preprocessingWDP(items,W);
		return 0;
	}
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 16);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		vector<GRBVar> newDim1 (nbItems);
		vector<GRBVar> newDim2 (nbItems);
		GRBVar newW = model.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS);
		
		// Initialization
		for (int i = 0; i < items.size(); i++){	
			newDim1[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			newDim2[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		}
		
		model.update();
		
		// Feasible/Infeasible patterns
		cout << "Print infeasible patterns" << endl;
		for (int i = 0; i < patternsF.size(); i++){
			GRBLinExpr sum = 0;
			for(int j = 0; j < patternsF[i].size();j++){
				if(patternsF[i][j] < nbItems)
					sum += newDim1[patternsF[i][j]];
				else 
					sum += newDim2[patternsF[i][j]-nbItems];
			}
			model.addConstr(sum <= newW);
			for(int j = 0; j < items.size();j++){
				bool isIP = false;
				for (int k = 0; k < patternsF[i].size(); k++){
					if (patternsF[i][k] == j || patternsF[i][k] == j + nbItems)
						isIP = true;		
				}
				if(!isIP && patternsC[i]+ items[j].lenO > W){
					/*cout << i << " " << patternsC[i] + items[j].widO << "-";
						for(int k = 0; k < patternsF[i].size();k++){
							cout << patternsF[i][k] << " ";
					}
					cout << j << endl;*/
					model.addConstr(sum + newDim1[j] >= newW + 1);
				}
				if(!isIP && patternsC[i]+ items[j].widO > W){
					/*cout << i << " " << patternsC[i] + items[j].lenO << "-";
						for(int k = 0; k < patternsF[i].size();k++){
							cout << patternsF[i][k] << " ";
					}
					cout << j + nbItems << endl;*/
					model.addConstr(sum + newDim2[j] >= newW + 1);			
				}
			}
		}
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();

		// Filling Solution
		cout << "W from " << W << " to " << newW.get(GRB_DoubleAttr_X) << endl;
		W = ceil(newW.get(GRB_DoubleAttr_X) - EPSILON);
		for (int i = 0; i < items.size(); i++){	
			cout << "Item " << i << " from " << items[i].lenO << " " << items[i].widO << " to " << newDim1[i].get(GRB_DoubleAttr_X) << " " << newDim2[i].get(GRB_DoubleAttr_X) << " ";
			items[i].widR = ceil(newDim1[i].get(GRB_DoubleAttr_X) - EPSILON);
			items[i].wid = ceil(newDim2[i].get(GRB_DoubleAttr_X) - EPSILON);
			cout << items[i].widR << " " << items[i].wid << endl;
		}	
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}


	// End
	return 0;
}

int preprocessingHLP(vector<Item>& items, int& H, const int& maxStack){
	int nbItems = items.size();
	// Local variables
	vector<vector<int> > patternsF; patternsF.push_back({});
	vector<int> patternsC; patternsC.push_back(0);
	for (int i = 0; i < items.size(); i++){	
		vector<vector<int> > patternsToAdd;
		for(int j = 0; j < patternsF.size();j++){
			if(patternsF[j].size() <= maxStack - 1 && patternsC[j] + items[i].heiO <= H){
				patternsToAdd.push_back(patternsF[j]); patternsToAdd.back().push_back(i);
				patternsC.push_back(patternsC[j] + items[i].heiO);
			}
		}
		for(int j = 0; j < patternsToAdd.size(); j++)
			patternsF.push_back(patternsToAdd[j]);
	}

	cout << "There are " << patternsF.size() << " feasible patterns " << endl;
	if(H <= 50){
		cout << "H is below 50, not much is to be gained by applying the LP preprocessing" << endl;
		preprocessingHDP(items,H);
		return 0;
	}
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 16);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		vector<GRBVar> newDim (nbItems);
		GRBVar newH = model.addVar(0, GRB_INFINITY, 1, GRB_CONTINUOUS);
		
		// Initialization
		for (int i = 0; i < items.size(); i++)
			newDim[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		
		model.update();
		
		// Feasible/Infeasible patterns
		for (int i = 0; i < patternsF.size(); i++){
			GRBLinExpr sum = 0;
			for(int j = 0; j < patternsF[i].size();j++)
				sum += newDim[patternsF[i][j]];
			model.addConstr(sum <= newH);
			if(patternsF[i].size() <= maxStack - 1){
				for(int j = 0; j < items.size();j++){
					bool isIP = false;
					for (int k = 0; k < patternsF[i].size(); k++){
						if (patternsF[i][k] == j)
							isIP = true;		
					}
					if(!isIP && patternsC[i]+ items[j].heiO > H)
						model.addConstr(sum + newDim[j] >= newH + 1);
				}
			}
		}
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();

		// Filling Solution
		cout << "H from " << H << " to " << newH.get(GRB_DoubleAttr_X) << endl;
		H = ceil(newH.get(GRB_DoubleAttr_X) - EPSILON);
		for (int i = 0; i < items.size(); i++){	
			cout << "Item " << i << " from " << items[i].heiO << " to " << newDim[i].get(GRB_DoubleAttr_X) << " ";
			items[i].hei = ceil(newDim[i].get(GRB_DoubleAttr_X) - EPSILON);
			cout << items[i].hei << endl;
		}	
	}

	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	// End
	return 0;
}

mycallbackKEP::mycallbackKEP(const Allocation& xallo, const vector<GRBVar>& xisItemPackedV, const vector<vector<GRBVar> >& xisArcUsed) {
    // initialize the Callback object
    this->ncuts = 0;  						// the number of cuts, initially 0
    this->allo = xallo;     				// the problem
	this->isItemPackedV = xisItemPackedV;   // the column starts
    this->isArcUsed = xisArcUsed;  			// the columns
}

int mycallbackKEP::getNcuts(){
	return this->ncuts; 
}

void mycallbackKEP::callback() {
    if (where == GRB_CB_MIPSOL) { 
        // identify the columns
		vector<int> succ (allo.nbItems,-1);
		vector<GRBLinExpr> succE (allo.nbItems,0);
		for (int i = 0; i < allo.cArcs.size(); i++){
			for (int j = 0; j < allo.maxStack-1; j++){
				if(ceil(getSolution(isArcUsed[i][j])- EPSILON) == 1){
					succ[allo.cArcs[i][0]] = allo.cArcs[i][1];
					succE[allo.cArcs[i][0]] += isArcUsed[i][j];
				}
			}
		}
		
		// for each column
		for (int i = 0; i < allo.nbItems; i++){
			if(ceil(getSolution(isItemPackedV[i])- EPSILON) == 1){
				cout << "Column " << i << " with items " << i << " ";
				GRBLinExpr LHS = 0; int RHS = 0; int height = allo.items[i].hei;
				int curr = i;
				while(succ[curr] != -1){ 
					LHS += succE[curr]; RHS += 1; height += allo.items[succ[curr]].hei; 
					curr = succ[curr];
					cout << curr << " ";
					if(height > allo.H){
						addLazy(LHS <= RHS -1); ncuts++; 			
						break;
					}
				}
				if (height > allo.H)
					cout << "height " << height << " -- forbidden!" << endl;
				else 
					cout << "height " << height << " -- allowed!" << endl;
			}
		}
	}
}

mycallbackYC::mycallbackYC(const Allocation& xallo, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoN, const vector<vector<vector<GRBVar> > >& xisItemCSPPackedInCoR, const vector<Item>& xitemsCSP, const vector<int>& xNPW) {
    // initialize the Callback object
    this->ncuts = 0;  										// the number of cuts, initially 0
    this->allo = xallo;     								// the problem
	this->isItemCSPPackedInCoN = xisItemCSPPackedInCoN;   	// the column starts
	this->isItemCSPPackedInCoR = xisItemCSPPackedInCoR;   	// the column starts	
    this->itemsCSP = xitemsCSP;  							// the columns
	this->NPW = xNPW;
	this->bestSol = -1;
}

int mycallbackYC::getNcuts(){
	return this->ncuts; 
}

vector<vector<vector<int> > > mycallbackYC::getOrdinates(){
	return this->ordinates; 
}

void mycallbackYC::callback() {
	if (where == GRB_CB_MIPSOL) {   
		cout << "*********** " << getDoubleInfo(GRB_CB_MIPSOL_OBJ) << "***************" << endl;
		if(floor(getDoubleInfo(GRB_CB_MIPSOL_OBJ) + EPSILON) <= bestSol){
			cout << "Already a solution like that" << endl;
			return;
		}
		
		int RHS = 0;
		GRBLinExpr LHS = 0;
		
		// Reduced instance (item type, W) 
		vector<int> isActivated(allo.W,0);
		vector<vector<int> > reducedInstance;
		for(int i = 0; i < isItemCSPPackedInCoN.size();i++){
			for(int j = 0; j < isItemCSPPackedInCoN[i].size();j++){
				for(int k = 0; k < isItemCSPPackedInCoN[i][j].size();k++){
					if(NPW[k] && k + itemsCSP[i].wid <= allo.W && ceil(getSolution(isItemCSPPackedInCoN[i][j][k])- EPSILON) == 1){
						LHS += isItemCSPPackedInCoN[i][j][k]; RHS++;
						for(int l = 0; l < pow(2,j);l++)
							reducedInstance.push_back({i,k,0});
						isActivated[k] = 1;
					}
					if(NPW[k] && k + itemsCSP[i].widR <= allo.W && ceil(getSolution(isItemCSPPackedInCoR[i][j][k])- EPSILON) == 1){
						LHS += isItemCSPPackedInCoR[i][j][k]; RHS++;
						for(int l = 0; l < pow(2,j);l++)
							reducedInstance.push_back({i,k,1});
						isActivated[k] = 1;
					}
				}
			}
		}
				
		int nbObj = reducedInstance.size();
		for(int i = 0; i < nbObj;i++){
			cout << reducedInstance[i][0] << " " << reducedInstance[i][1] << " " << reducedInstance[i][2] << endl;
		}
	
		// Local variables
		IloEnv env;		
		
		// Model
		try{
			IloModel model(env);
			
			IloIntervalVarArray Tasks(env);
			IloIntervalVarArray2 BinTasks(env, allo.W);
			for (int i = 0; i < allo.W; i++)
				BinTasks[i] = IloIntervalVarArray(env);

			// Init Tasks
			for (int i = 0; i < nbObj; i++) {
				if (reducedInstance[i][2] == 0){
					IloIntervalVar ADD = IloIntervalVar(env, itemsCSP[reducedInstance[i][0]].len);
					int sm = 0; int sM = allo.L - itemsCSP[reducedInstance[i][0]].len; int em = itemsCSP[reducedInstance[i][0]].len; int eM = allo.L;
					for (int j = 0; j < nbObj; j++) {
						if (reducedInstance[i] == reducedInstance[j]){
							if (i > j){
								sm += itemsCSP[reducedInstance[i][0]].len;
								em += itemsCSP[reducedInstance[i][0]].len;
							}
							if (i < j){
								sM -= itemsCSP[reducedInstance[i][0]].len;
								eM -= itemsCSP[reducedInstance[i][0]].len;
							}
						}
					}
					ADD.setEndMax(eM); ADD.setEndMin(em);
					ADD.setStartMin(sm); ADD.setStartMax(sM);
					Tasks.add(ADD);
					for (int j = 0; j < itemsCSP[reducedInstance[i][0]].wid; j++){
						if(isActivated[reducedInstance[i][1] + j])
							BinTasks[reducedInstance[i][1] + j].add(Tasks[i]);
					}
				}
				else {
					IloIntervalVar ADD = IloIntervalVar(env, itemsCSP[reducedInstance[i][0]].lenR);
					int sm = 0; int sM = allo.L - itemsCSP[reducedInstance[i][0]].lenR; int em = itemsCSP[reducedInstance[i][0]].lenR; int eM = allo.L;
					for (int j = 0; j < nbObj; j++) {
						if (reducedInstance[i] == reducedInstance[j]){
							if (i > j){
								sm += itemsCSP[reducedInstance[i][0]].lenR;
								em += itemsCSP[reducedInstance[i][0]].lenR;
							}
							if (i < j){
								sM -= itemsCSP[reducedInstance[i][0]].lenR;
								eM -= itemsCSP[reducedInstance[i][0]].lenR;
							}
						}
					}
					ADD.setEndMax(eM); ADD.setEndMin(em);
					ADD.setStartMin(sm); ADD.setStartMax(sM);
					Tasks.add(ADD);
					for (int j = 0; j < itemsCSP[reducedInstance[i][0]].widR; j++){
						if(isActivated[reducedInstance[i][1] + j])
							BinTasks[reducedInstance[i][1] + j].add(Tasks[i]);
					}
				}
			}

			// Remove symmetry
			for (int i = 0; i < nbObj; i++) {
				for (int j = 0; j < i; j++) {
					if (reducedInstance[i] == reducedInstance[j])
						model.add(IloEndBeforeStart(env, Tasks[j], Tasks[i]));
				}
			}

			for(int i = 0; i < allo.W;i++){
				if(isActivated[i])
					model.add(IloNoOverlap(env, BinTasks[i]));
			}		

			IloCP cp(model);
			cp.setParameter(IloCP::NoOverlapInferenceLevel, IloCP::Extended);
			cp.setParameter(IloCP::IntervalSequenceInferenceLevel, IloCP::Extended);
			cp.setOut(env.getNullStream());

			// If infeasible
			bool check = cp.solve();
			if (!check){
				if (cp.getStatus() == 3){ // If infeasible
					cout << "Infeasible" << endl;
					addLazy(LHS <= RHS -1); ncuts++; 	
					return;
				}
			}
			else{
				bestSol = floor(getDoubleInfo(GRB_CB_MIPSOL_OBJ) + EPSILON);
				ordinates.resize(0); ordinates.resize(itemsCSP.size());	
				for (int i = 0; i < nbObj; i++)
					ordinates[reducedInstance[i][0]].push_back({cp.getStart(Tasks[i]),reducedInstance[i][1],reducedInstance[i][2]});
			}
		}

		// Exceptions
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Exception during optimization" << endl;
		}
	}
}		
