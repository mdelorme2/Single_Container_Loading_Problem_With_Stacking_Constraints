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
	model(allo);
	allo.rescaleW();
	allo.printInfo(pathAndFileout);
	allo.printFile(picture);
}

int model(Allocation& allo){

	// Transform into CSP
	vector<Item> itemsCSP;
	vector<int> demandCSP;
	vector<vector<int> > indexCSP;
	for(int i =0; i<allo.items.size();i++){
		bool exist = false;
		for(int j = 0; j<itemsCSP.size();j++){
			if(allo.items[i].len == itemsCSP[j].len && allo.items[i].wid == itemsCSP[j].wid && allo.items[i].lenR == itemsCSP[j].lenR && allo.items[i].widR == itemsCSP[j].widR){
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
	
	// Local variables
	vector<int> NPL (allo.L+1,0);
	vector<int> NPW (allo.W+1,0);
	int minW = allo.W;
	int minL = allo.L;
	
	// Compute normal patterns
	NPL[0] = 1; NPW[0] = 1;
	for(int i =0; i<indexCSP.size();i++){
		minW = min(minW,itemsCSP[i].wid); minW = min(minW,itemsCSP[i].widR); 
		minL = min(minL,itemsCSP[i].len); minL = min(minL,itemsCSP[i].lenR); 
		for(int j = 0; j < demandCSP[i]; j++){
			for(int k = allo.L - min(itemsCSP[i].len,itemsCSP[i].lenR); k>=0; k--){
				if(NPL[k] == 1 && k + itemsCSP[i].len <= allo.L) NPL[k + itemsCSP[i].len] = 1;
				if(NPL[k] == 1 && k + itemsCSP[i].lenR <= allo.L) NPL[k + itemsCSP[i].lenR] = 1;
			}
			for(int k = allo.W - min(itemsCSP[i].wid,itemsCSP[i].widR); k>=0; k--){
				if(NPW[k] == 1 && k + itemsCSP[i].wid <= allo.W) NPW[k + itemsCSP[i].wid] = 1;
				if(NPW[k] == 1 && k + itemsCSP[i].widR <= allo.W) NPW[k + itemsCSP[i].widR] = 1;
			}
		}
	}
	
	// Remove ending NPs
	for(int i = allo.W - minW + 1; i < allo.W; i++) NPW[i] = 0;
	for(int i = allo.L - minL + 1; i < allo.L; i++) NPL[i] = 0;
	
	// Print normal patterns
	cout << "MinL and minW are " << minL << " " << minW << endl;
	cout << "Normal patterns activated " << endl;
	cout << "L:";
	for(int i = 0; i <= allo.L; i++){
		if (NPL[i]) cout << i << " ";
	}
	cout << endl << "W:";
	for(int i = 0; i <= allo.W; i++){
		if (NPW[i]) cout << i << " ";
	}	
	cout << endl;	
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 50);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr totWeight = 0;
		vector<vector<GRBLinExpr> > isSpaceBusy(allo.L);
		
		// Related to the item type 
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoN (itemsCSP.size());
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoR (itemsCSP.size());
		vector<GRBVar> isItemCSPPackedV(itemsCSP.size());		
		vector<GRBLinExpr> isItemCSPPacked(itemsCSP.size(),0);
		vector<GRBLinExpr> isItemPacked(itemsCSP.size(),0);
		
		// Related to the item 
		vector<vector<vector<GRBVar> > > isItemStackedInCo (allo.nbItems);
		vector<GRBLinExpr> isItemStacked(allo.nbItems,0);;
		vector<GRBLinExpr> heighColumn(allo.nbItems,0);
		vector<GRBVar> isItemStackedV(allo.nbItems);
		vector<GRBVar> isItemPackedV(allo.nbItems);
		
		// Initialization
		for (int i = 0; i < itemsCSP.size(); i++){	
			isItemCSPPackedV[i] = model.addVar(0, demandCSP[i], 0, GRB_INTEGER);
			isItemCSPPackedInCoN[i].resize(allo.L);
			isItemCSPPackedInCoR[i].resize(allo.L);
			for (int j = 0; j < allo.L; j++){	
				if(NPL[j] && j + itemsCSP[i].len <= allo.L){
					isItemCSPPackedInCoN[i][j].resize(allo.W);
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].wid <= allo.W){
							if((j + itemsCSP[i].len != allo.L && j + itemsCSP[i].len + minL > allo.L && NPL[allo.L - itemsCSP[i].len]) || (k + itemsCSP[i].wid != allo.W && k + itemsCSP[i].wid + minW > allo.W && NPW[allo.W - itemsCSP[i].wid])){
								isItemCSPPackedInCoN[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
							}
							else{
								isItemCSPPackedInCoN[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);
							}
						}
					}
				}
				if(NPL[j] && j + itemsCSP[i].lenR <= allo.L){
					isItemCSPPackedInCoR[i][j].resize(allo.W);
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].widR <= allo.W){
							if(itemsCSP[i].lenR == itemsCSP[i].len && itemsCSP[i].widR == itemsCSP[i].wid)
								isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
							else{
								if((j + itemsCSP[i].lenR != allo.L  && j + itemsCSP[i].lenR + minL > allo.L && NPL[allo.L - itemsCSP[i].lenR]) || (k + itemsCSP[i].widR != allo.W && k + itemsCSP[i].widR + minW > allo.W && NPW[allo.W - itemsCSP[i].widR]))
									isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 0, 0, GRB_BINARY);
								else
									isItemCSPPackedInCoR[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);
							}								
						}
					}
				}
			}
		}
		for (int i = 0; i < allo.items.size(); i++){	
			isItemStackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemStackedInCo[i].resize(allo.items.size());
			for (int j = 0; j < allo.items.size(); j++){
				if(allo.cAdj[i][j]){
					isItemStackedInCo[i][j].resize(allo.maxStack);
					for (int k = 1; k < allo.maxStack; k++)
						isItemStackedInCo[i][j][k] = model.addVar(0, 1, 0, GRB_BINARY);
				}
			}
		}

		for (int j = 0; j < allo.L; j++)
			isSpaceBusy[j].resize(allo.W,0);
			
		model.update();
		
		// Compute values packing
		for (int i = 0; i < itemsCSP.size(); i++){	
			for (int j = 0; j < allo.L; j++){	
				if(NPL[j] && j + itemsCSP[i].len <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].wid <= allo.W){
							isItemCSPPacked[i] += isItemCSPPackedInCoN[i][j][k];
							for(int j2 = 0; j2 < itemsCSP[i].len ; j2++){
								if(NPL[j+j2]){
									for(int k2 = 0; k2 < itemsCSP[i].wid ; k2++){
										if(NPW[k+k2])
											isSpaceBusy[j+j2][k+k2] += isItemCSPPackedInCoN[i][j][k];
									}
								}
							}
						}
					}
				}
				if(NPL[j] && j + itemsCSP[i].lenR <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].widR <= allo.W){
							isItemCSPPacked[i] += isItemCSPPackedInCoR[i][j][k];
							for(int j2 = 0; j2 < itemsCSP[i].lenR ; j2++){
								if(NPL[j+j2]){
									for(int k2 = 0; k2 < itemsCSP[i].widR ; k2++){
										if(NPW[k+k2])
											isSpaceBusy[j+j2][k+k2] += isItemCSPPackedInCoR[i][j][k];
									}
								}
							}
						}
					}
				}
			}
			for(int j = 0; j < indexCSP[i].size();j++)
				isItemPacked[i] += isItemPackedV[indexCSP[i][j]];
			model.addConstr(isItemCSPPacked[i] == isItemCSPPackedV[i]);
		}
	
		// Constraints -- stacking
		for (int i = 0; i < allo.items.size(); i++){
			for (int k = 1; k < allo.maxStack; k++){
				GRBLinExpr isPosBusy = 0;
				// Stack position 2
				if(k == 1){
					for (int j = 0; j < allo.items.size(); j++){
						if(allo.cAdj[i][j]){
							model.addConstr(isItemStackedInCo[i][j][1] <= isItemPackedV[i]);
							isPosBusy += isItemStackedInCo[i][j][1];
							isItemStacked[j] += isItemStackedInCo[i][j][1];
							heighColumn[i] += isItemStackedInCo[i][j][1] * allo.items[j].hei;
						}
					}
				}
				else{
					// Stack position 3 or above
					for (int j = 0; j < allo.items.size(); j++){
						if(allo.cAdj[i][j]){
							GRBLinExpr tSum = 0;
							for (int l = 0; l < allo.items.size(); l++){
								if(allo.cAdj[i][l] && allo.cAdj[l][j]){
									tSum += isItemStackedInCo[i][l][k-1];
								}
							}
							model.addConstr(isItemStackedInCo[i][j][k] <= tSum);
							isPosBusy += isItemStackedInCo[i][j][k];
							isItemStacked[j] += isItemStackedInCo[i][j][k];
							heighColumn[i] += isItemStackedInCo[i][j][k] * allo.items[j].hei;
						}				
					}
				}
				model.addConstr(isPosBusy <= 1);
			}
		}

		// Global constraints -- linking var between stacking and packing
		for (int i = 0; i < itemsCSP.size(); i++)
			model.addConstr(isItemPacked[i] == isItemCSPPackedV[i]);
	
		// Global constraints -- linking expr and var for stacking
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(isItemStacked[i] == isItemStackedV[i]);
	
		// Global constraints -- weight
		for (int i = 0; i < allo.items.size(); i++)
			totWeight += (isItemPackedV[i] + isItemStackedV[i]) * allo.items[i].wei;
		model.addConstr(totWeight <= allo.maxWeight);

		// Global constraints -- height
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(heighColumn[i] <= (allo.H - allo.items[i].hei) * isItemPackedV[i]);				
		
		// Global constraints -- one item at most
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(isItemPackedV[i] + isItemStackedV[i] <= 1);

		// Objective function -- maximize volume used
		for (int i = 0; i < allo.items.size(); i++)
			objFun += (isItemPackedV[i] + isItemStackedV[i]) * allo.items[i].pro;		
		
		model.setObjective(objFun, GRB_MAXIMIZE);
		
		// Nooverlap
		for (int k = 0; k < allo.L; k++){	
			if(NPL[k]){
				for (int l = 0; l < allo.W; l++){
					if(NPW[l])
						model.addConstr(isSpaceBusy[k][l] <= 1);
				}
			}
		}
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
		model.optimize();
		
		// Filling Info
		allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
		allo.infos.UB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
		allo.infos.opt = false;
		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);

		// If no solution found
		if (model.get(GRB_IntAttr_SolCount) < 1){
			cout << "Failed to optimize ILP. " << endl;
			allo.infos.LB  = 0;
			return -1;
		}

		// If solution found
		allo.infos.LB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);	
		if(allo.infos.LB == allo.infos.UB) allo.infos.opt = true;

		// Filling Solution
		allo.pps.resize(allo.items.size());
		allo.columns.resize(allo.items.size());
		for (int i = 0; i < itemsCSP.size(); i++){	
			cout << i << " " << ceil(isItemCSPPackedV[i].get(GRB_DoubleAttr_X) - EPSILON) << endl;
			for (int j = 0; j < allo.L; j++){	
				if(NPL[j] && j + itemsCSP[i].len <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].wid <= allo.W && ceil(isItemCSPPackedInCoN[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
 							int index = indexCSP[i].back();
							while(ceil(isItemPackedV[index].get(GRB_DoubleAttr_X) - EPSILON) == 0){
								indexCSP[i].pop_back();
								index = indexCSP[i].back();
							}
							allo.pps[index].push_back(j); allo.pps[index].push_back(k); allo.pps[index].push_back(0); 
							allo.columns[index].push_back(index); 
							for(int m = 1; m<allo.maxStack;m++){
								for (int l = 0; l < allo.items.size(); l++){
									if(allo.cAdj[index][l] && ceil(isItemStackedInCo[index][l][m].get(GRB_DoubleAttr_X) - EPSILON) == 1)
										allo.columns[index].push_back(l); 
								}
							}
							indexCSP[i].pop_back();
						}
					}
				}
			}
			for (int j = 0; j < allo.L; j++){	
				if(NPL[j] && j + itemsCSP[i].lenR <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(NPW[k] && k + itemsCSP[i].widR <= allo.W && ceil(isItemCSPPackedInCoR[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
							int index = indexCSP[i].back();
							while(ceil(isItemPackedV[index].get(GRB_DoubleAttr_X) - EPSILON) == 0){
								indexCSP[i].pop_back();
								index = indexCSP[i].back();
							}
							allo.pps[index].push_back(j); allo.pps[index].push_back(k); allo.pps[index].push_back(1); 
							allo.columns[index].push_back(index); 
							for(int m = 1; m<allo.maxStack;m++){
								for (int l = 0; l < allo.items.size(); l++){
									if(allo.cAdj[index][l] && ceil(isItemStackedInCo[index][l][m].get(GRB_DoubleAttr_X) - EPSILON) == 1)
										allo.columns[index].push_back(l); 
								}
							}
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
	env.set(GRB_DoubleParam_MemLimit, 50);
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
	env.set(GRB_DoubleParam_MemLimit, 50);
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
