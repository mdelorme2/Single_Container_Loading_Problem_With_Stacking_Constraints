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
	
	// main loop
	allo.infos.nbCuts2 = 0;
	vector<vector<vector<int> > > cuts1; 
	vector<vector<vector<int> > > cuts2;
	vector<Item> itemsCSP; 
	vector<int> demandCSP;
	vector<vector<int> > indexCSP;
	int ret1, ret2, ret3;
	step0(allo,itemsCSP,demandCSP,indexCSP);
	allo.infos.nbCuts1 = 0;
	allo.infos.nbCuts2 = 0;
	for(int i = 0; i < 3; i++){
		if(i == 0){
			cuts2.resize(0);
			ret1 = step1(allo,cuts1,itemsCSP,demandCSP,indexCSP);
			if (ret1 == 2) break;
			continue;
		}
		if(i == 1){
			ret2 = step2(allo,cuts1.back(),cuts2,itemsCSP,indexCSP);
			if (ret2 == 2) break;
			if(ret2 == 1){
				i = -1; 
				allo.infos.nbCuts1++;
				continue;
			}
		}
		if(i == 2){
			ret3 = step3(allo,cuts2.back(),itemsCSP,indexCSP);
			if(ret3 == 1){
				i = 0; 
				allo.infos.nbCuts2++;
				continue;
			}
		}
	}
	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
	if(ret1 == 0 && ret2 == 0 && ret3 == 0){
		allo.rescaleW();
		allo.printInfo(pathAndFileout);
		allo.printFile(picture);
	}
	else{
		allo.infos.opt = 0;
		allo.infos.LB = 0;
		allo.printInfo(pathAndFileout);
	}
}

int step0(Allocation& allo, vector<Item>& itemsCSP, vector<int>& demandCSP, vector<vector<int> >& indexCSP){

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
	return 0;
}

int step1(Allocation& allo, vector<vector<vector<int> > > & cuts1, const vector<Item>& itemsCSP, const vector<int>& demandCSP, const vector<vector<int> >& indexCSP){

	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 50);
	env.start();
	
	// Step 1 -- make columns
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr totWeight = 0;
		GRBLinExpr totArea = 0;

		// Related to the item type 
		vector<vector<GRBVar> > isItemCSPPackedV(itemsCSP.size());		
		vector<GRBLinExpr> isItemPacked(itemsCSP.size(),0);
		
		// Related to the item 
		vector<GRBVar> isItemPackedV(allo.nbItems);
		vector<GRBVar> isItemStackedV(allo.nbItems);
		vector<vector<vector<GRBVar> > > isItemStackedInCo (allo.nbItems);
		vector<GRBLinExpr> isItemStacked(allo.nbItems,0);
		vector<GRBLinExpr> heighColumn(allo.nbItems,0);

		// Initialization
		for (int i = 0; i < itemsCSP.size(); i++){
			isItemCSPPackedV[i].resize(ceil(log2(demandCSP[i]))+1);
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				isItemCSPPackedV[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
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
		
		model.update();

		// Compute values packing
		for (int i = 0; i < itemsCSP.size(); i++){			
			for(int j = 0; j < indexCSP[i].size();j++)
				isItemPacked[i] += isItemPackedV[indexCSP[i][j]];
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

		// Global constraints -- weight 
		for (int i = 0; i < allo.items.size(); i++)
			totWeight += (isItemPackedV[i] + isItemStackedV[i]) * allo.items[i].wei;
		model.addConstr(totWeight <= allo.maxWeight);

		// Global constraints -- area
		for (int i = 0; i < itemsCSP.size(); i++){			
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++)
				totArea += pow(2,j) * isItemCSPPackedV[i][j] * min(itemsCSP[i].len * itemsCSP[i].wid, itemsCSP[i].lenR * itemsCSP[i].widR);	
		}
		model.addConstr(totArea <= allo.W * allo.L);
		
		// Global constraints -- linking var between stacking and packing
		for (int i = 0; i < itemsCSP.size(); i++){
			GRBLinExpr number = 0;
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++) 
				number += pow(2,j) * isItemCSPPackedV[i][j];
			model.addConstr(isItemPacked[i] == number);
		}

		// Global constraints -- linking expr and var
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(isItemStacked[i] == isItemStackedV[i]);
		
		// Global constraints -- height
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(heighColumn[i] <= (allo.H - allo.items[i].hei) * isItemPackedV[i]);			
		
		// Global constraints -- one item at most
		for (int i = 0; i < allo.items.size(); i++)
			model.addConstr(isItemPackedV[i] + isItemStackedV[i] <= 1);

		// Global constraints -- no good cuts
		for (int i = 0; i < cuts1.size(); i++){
			int RHS = 0;
			GRBLinExpr LHS = 0;
			for(int j = 0; j < cuts1[i].size(); j++){
				RHS++;
				LHS += isItemCSPPackedV[cuts1[i][j][0]][cuts1[i][j][1]];
			}
			model.addConstr(LHS <= RHS - 1);
		}
			
		// Objective function -- maximize volume used
		for (int i = 0; i < allo.items.size(); i++)
			objFun += (isItemPackedV[i] + isItemStackedV[i]) * allo.items[i].pro;		

		model.setObjective(objFun, GRB_MAXIMIZE);
	
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  max(0.0,3600 - (getCPUTime() - initTimeModelCPU)));
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
		model.optimize();
		
		// Filling Info
		allo.infos.LB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		allo.infos.UB = floor(model.get(GRB_DoubleAttr_ObjBound) + EPSILON);
		if(allo.infos.LB == allo.infos.UB) allo.infos.opt = true;
		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);

		if(model.get(GRB_IntAttr_Status) == 9){
			return 2;
		}
		
		// Filling Info and solution
		allo.columns.resize(0);
		allo.columns.resize(allo.items.size());
		for (int i = 0; i < allo.items.size(); i++){	
			if(ceil(isItemPackedV[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				allo.columns[i].push_back(i);		
				for(int m = 1; m<allo.maxStack;m++){
					for (int l = 0; l < allo.items.size(); l++){
						if(allo.cAdj[i][l] && ceil(isItemStackedInCo[i][l][m].get(GRB_DoubleAttr_X) - EPSILON) == 1)
							allo.columns[i].push_back(l); 
					}
				}
			}
		}	

		cout << "Solution step 1:[";
		vector<vector<int> > newCut;
		for (int i = 0; i < itemsCSP.size(); i++){			
			int number = 0;
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				if(ceil(isItemCSPPackedV[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
					number += pow(2,j);
					newCut.push_back({i,j});
				}
			}
			cout << number << " ";
		}
		cout << "]" << endl;
		cuts1.push_back(newCut);
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

int step2(Allocation& allo, const vector<vector<int> > & cut1, vector<vector<vector<int> > > & cuts2, const vector<Item>& itemsCSP, const vector<vector<int> >& indexCSP){

	// Demand CSP 
	vector<int> demandCSP (itemsCSP.size(),0);
	for(int i = 0; i < cut1.size();i++)
		demandCSP[cut1[i][0]] += pow(2,cut1[i][1]);
		
	// Local variables
	int W = allo.W;
	vector<int> NPW (W+1,0);
	int minW = W;
	int nbItems = itemsCSP.size();
	
	// Compute normal patterns
	NPW[0] = 1;
	for(int i =0; i<nbItems;i++){
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
	
	// Print normal patterns
	/*cout << "minW is " << minW << endl;
	cout << "Normal patterns for W activated " << endl;
	for(int i = 0; i <= allo.W; i++){
		if (NPW[i]) cout << i << " ";
	}	
	cout << endl;	*/
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 50);
	env.start();

	// Step 1 -- make columns
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		vector<GRBLinExpr> columnFilled (W,0);

		// Related to the item type 
		vector<GRBVar> isItemCSPPackedV(itemsCSP.size());			
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoN (itemsCSP.size());
		vector<vector<vector<GRBVar> > > isItemCSPPackedInCoR (itemsCSP.size());
		vector<GRBLinExpr> isItemCSPPacked(itemsCSP.size(),0);		
		vector<GRBLinExpr> isItemPacked(itemsCSP.size(),0);

		// Initialization
		for (int i = 0; i < itemsCSP.size(); i++){
			if(demandCSP[i] == 0) continue;
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
			model.addConstr(isItemCSPPacked[i] == isItemCSPPackedV[i]);
		}

		// Global constraints -- length
		for (int i = 0; i < allo.W; i++){
			if(NPW[i])
				model.addConstr(columnFilled[i] <= allo.L);
		}
		
		// Global constraints -- satisfy demand
		for (int i = 0; i < itemsCSP.size(); i++){
			if(demandCSP[i] == 0) continue;
			model.addConstr(isItemCSPPackedV[i] == demandCSP[i]);
		}

		// Global constraints -- no good cuts
		for (int i = 0; i < cuts2.size(); i++){
			int RHS = 0;
			GRBLinExpr LHS = 0;
			for(int j = 0; j < cuts2[i].size(); j++){
				RHS++;
				if(cuts2[i][j][0] == 0) LHS += isItemCSPPackedInCoN[cuts2[i][j][1]][cuts2[i][j][2]][cuts2[i][j][3]];
				if(cuts2[i][j][0] == 1) LHS += isItemCSPPackedInCoR[cuts2[i][j][1]][cuts2[i][j][2]][cuts2[i][j][3]];
			}
			model.addConstr(LHS <= RHS - 1);
		}
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  max(1.0,3600 - (getCPUTime() - initTimeModelCPU)));
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
	//	model.getEnv().set(GRB_IntParam_OutputFlag,0);
		model.optimize();
		
		// If infeasible
		if(model.get(GRB_IntAttr_Status) == 3)
			return 1;

		allo.infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		allo.infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		allo.infos.nbNZ = model.get(GRB_IntAttr_NumNZs);	

		// If time limit
		if(model.get(GRB_IntAttr_Status) == 9)
			return 2;
		
		vector<vector<int> > newCut;
		for (int i = 0; i < itemsCSP.size(); i++){		
			if(demandCSP[i] == 0) continue;		
			for (int j = 0; j < ceil(log2(demandCSP[i]))+1; j++){
				for (int k = 0; k < allo.W; k++){
					if(NPW[k] && k + itemsCSP[i].wid <= allo.W){
						if(ceil(isItemCSPPackedInCoN[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1)
							newCut.push_back({0,i,j,k});
					}
					if(NPW[k] && k + itemsCSP[i].widR <= allo.W){
						if(ceil(isItemCSPPackedInCoR[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1)
							newCut.push_back({1,i,j,k});
					}
				}
			}
		}
		cuts2.push_back(newCut);
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
	
int step3(Allocation& allo, const vector<vector<int> > & cut2, const vector<Item>& itemsCSP, const vector<vector<int> >& indexCSP){

	// Reduced instance (item type, W) 
	vector<int> isActivated(allo.W,0);
	vector<vector<int> > reducedInstance;
	for(int i = 0; i < cut2.size();i++){
		isActivated[cut2[i][3]] = 1;
		for(int j = 0; j < pow(2,cut2[i][2]);j++)
			reducedInstance.push_back({cut2[i][1],cut2[i][3], cut2[i][0]});
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
		cp.setParameter(IloCP::TimeLimit, max(1.0,3600 - (getCPUTime() - initTimeModelCPU)));
	//	cp.setOut(env.getNullStream());

		// If infeasible
		bool check = cp.solve();
		if (!check){
			env.error() << "Failed to optimize CP. " << cp.getStatus() << endl;
			if (cp.getStatus() == 3){ // If infeasible
				env.end();
				return 1;
			}
			else { // If time limit
				cout << "Too long" << endl;
				env.end();
				return 2;
			}
		}
		
		// Filling Solution
		allo.pps.resize(allo.items.size());
		vector<int> isItemNotDone (allo.items.size(),0);
		for(int i = 0; i < allo.items.size();i++){
			if(allo.columns[i].size()>0) 
				isItemNotDone[i] = 1;
		}
		for (int i = 0; i < nbObj; i++){	
			// Find item index
			int index = -1;
			for(int j = 0; j < indexCSP[reducedInstance[i][0]].size();j++){
				if(isItemNotDone[indexCSP[reducedInstance[i][0]][j]] == 1){
					index = indexCSP[reducedInstance[i][0]][j];
					isItemNotDone[index] = 0;
					cout << "pack item " << i << " which is in reality item " << index << endl;
					break;
				}
			}
			if(reducedInstance[i][2] == 0){
				allo.pps[index].push_back(cp.getStart(Tasks[i])); allo.pps[index].push_back(reducedInstance[i][1]); allo.pps[index].push_back(0); 
			}
			if(reducedInstance[i][2] == 1){
				allo.pps[index].push_back(cp.getStart(Tasks[i])); allo.pps[index].push_back(reducedInstance[i][1]); allo.pps[index].push_back(1);
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
	//cout << "W from " << W << " to ";
	while (!isR[W]) W--;
	//cout << W << endl;
	
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
		//cout << "Item " << i << " wid from " << items[i].wid << " to ";
		while(!isR[W - items[i].wid]) items[i].wid++;
		//cout << items[i].wid << endl;
		//cout << "Item " << i << " widR from " << items[i].widR << " to ";
		while(!isR[W - items[i].widR]) items[i].widR++;
		//cout << items[i].widR << endl;
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
	//cout << "L from " << L << " to ";
	while (!isR[L]) L--;
	//cout << L << endl;
	
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
		//cout << "Item " << i << " len from " << items[i].len << " to ";
		while(!isR[L - items[i].len]) items[i].len++;
		//cout << items[i].len << endl;
		//cout << "Item " << i << " lenR from " << items[i].lenR << " to ";
		while(!isR[L - items[i].lenR]) items[i].lenR++;
		//cout << items[i].lenR << endl;
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
	//cout << "H from " << H << " to ";
	while (!isR[H]) H--;
	//cout << H << endl;
	
	// Preprocessing 2 -- increase height
	for (int i = 0; i < nbItems; i++){
		isR.resize(0); isR.resize(H+1,false); isR[0] = true;
		for (int j = 0; j < nbItems; j++){
			if(j != i){
				for(int k = H - items[j].hei; k>=0; k--)
					if(isR[k] == 1) isR[k + items[j].hei] = 1;
			}
		}
		//cout << "Item " << i << " height from " << items[i].hei << " to ";
		while(!isR[H - items[i].hei]) items[i].hei++;
		//cout << items[i].hei << endl;
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
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.optimize();

		// Filling Solution
		cout << "W from " << W << " to " << newW.get(GRB_DoubleAttr_X) << endl;
		W = ceil(newW.get(GRB_DoubleAttr_X) - EPSILON);
		for (int i = 0; i < items.size(); i++){	
			//cout << "Item " << i << " from " << items[i].lenO << " " << items[i].widO << " to " << newDim1[i].get(GRB_DoubleAttr_X) << " " << newDim2[i].get(GRB_DoubleAttr_X) << " ";
			items[i].widR = ceil(newDim1[i].get(GRB_DoubleAttr_X) - EPSILON);
			items[i].wid = ceil(newDim2[i].get(GRB_DoubleAttr_X) - EPSILON);
			//cout << items[i].widR << " " << items[i].wid << endl;
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
		model.getEnv().set(GRB_IntParam_OutputFlag, 0);
		model.optimize();

		// Filling Solution
		//cout << "H from " << H << " to " << newH.get(GRB_DoubleAttr_X) << endl;
		H = ceil(newH.get(GRB_DoubleAttr_X) - EPSILON);
		for (int i = 0; i < items.size(); i++){	
			// cout << "Item " << i << " from " << items[i].heiO << " to " << newDim[i].get(GRB_DoubleAttr_X) << " ";
			items[i].hei = ceil(newDim[i].get(GRB_DoubleAttr_X) - EPSILON);
			// cout << items[i].hei << endl;
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

