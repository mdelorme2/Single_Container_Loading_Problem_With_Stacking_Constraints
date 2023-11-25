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
	model(allo);
	allo.printInfo(pathAndFileout);
	allo.printFile(picture);
}

int model(Allocation& allo){

	allo.infos.timeCPU.push_back(getCPUTime() - initTimeModelCPU);
	
	GRBEnv env = GRBemptyenv;
	env.set(GRB_DoubleParam_MemLimit, 50);
	env.start();
		
	// Model
	try{
		// Local variables
		GRBModel model = GRBModel(env);
		GRBLinExpr objFun = 0;
		GRBLinExpr totWeight = 0;
		vector<GRBVar> isItemPackedInCoN (allo.nbItems);
		vector<GRBVar> isItemPackedInCoR (allo.nbItems);
		vector<vector<vector<GRBVar> > > isItemStackedInCo (allo.nbItems);
		vector<GRBLinExpr> isItemStacked(allo.nbItems,0);
		vector<GRBLinExpr> isItemPacked(allo.nbItems,0);
		vector<GRBLinExpr> heighColumn(allo.nbItems,0);
		vector<GRBVar> isItemStackedV(allo.nbItems);
		vector<GRBVar> isItemPackedV(allo.nbItems);
		vector<GRBVar> isItemPackedInCoX (2*allo.nbItems);
		vector<GRBVar> isItemPackedInCoY (2*allo.nbItems);
		vector<vector<GRBVar> > aij (2*allo.nbItems);
		vector<vector<GRBVar> > bij (2*allo.nbItems);
		vector<vector<GRBVar> > cij (2*allo.nbItems);
		vector<vector<GRBVar> > dij (2*allo.nbItems);

		// Initialization
		for (int i = 0; i < allo.items.size(); i++){	
			isItemStackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedInCoN[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedInCoR[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemStackedInCo[i].resize(allo.items.size());
			isItemPackedInCoX[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS); isItemPackedInCoX[i + allo.nbItems] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			isItemPackedInCoY[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS); isItemPackedInCoY[i + allo.nbItems] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			aij[i].resize(2*allo.nbItems); aij[i+allo.nbItems].resize(2*allo.nbItems);
			bij[i].resize(2*allo.nbItems); bij[i+allo.nbItems].resize(2*allo.nbItems);
			cij[i].resize(2*allo.nbItems); cij[i+allo.nbItems].resize(2*allo.nbItems);
			dij[i].resize(2*allo.nbItems); dij[i+allo.nbItems].resize(2*allo.nbItems);			
			for (int j = 0; j < allo.items.size(); j++){
				if(allo.cAdj[i][j]){
					isItemStackedInCo[i][j].resize(allo.maxStack);
					for (int k = 1; k < allo.maxStack; k++)
						isItemStackedInCo[i][j][k] = model.addVar(0, 1, 0, GRB_BINARY);
				}
				aij[i][j] = model.addVar(0, 1, 0, GRB_BINARY); aij[i][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY); aij[i+allo.nbItems][j] = model.addVar(0, 1, 0, GRB_BINARY); aij[i+allo.nbItems][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY);
				bij[i][j] = model.addVar(0, 1, 0, GRB_BINARY); bij[i][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY); bij[i+allo.nbItems][j] = model.addVar(0, 1, 0, GRB_BINARY); bij[i+allo.nbItems][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY);
				cij[i][j] = model.addVar(0, 1, 0, GRB_BINARY); cij[i][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY); cij[i+allo.nbItems][j] = model.addVar(0, 1, 0, GRB_BINARY); cij[i+allo.nbItems][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY);
				dij[i][j] = model.addVar(0, 1, 0, GRB_BINARY); dij[i][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY); dij[i+allo.nbItems][j] = model.addVar(0, 1, 0, GRB_BINARY); dij[i+allo.nbItems][j+allo.nbItems] = model.addVar(0, 1, 0, GRB_BINARY);	
			}
		}

		model.update();
		
		// Compute values packing
		for (int i = 0; i < allo.items.size(); i++){
			model.addConstr(isItemPackedInCoX[i] + allo.items[i].wid <= allo.W + allo.W * (1 - isItemPackedInCoN[i])); 
			model.addConstr(isItemPackedInCoX[i + allo.nbItems] + allo.items[i].widR <= allo.W + allo.W * (1 - isItemPackedInCoR[i])); 
			model.addConstr(isItemPackedInCoY[i] + allo.items[i].len <= allo.L + allo.L * (1 - isItemPackedInCoN[i])); 
			model.addConstr(isItemPackedInCoY[i + allo.nbItems] + allo.items[i].lenR <= allo.L + allo.L * (1 - isItemPackedInCoR[i])); 
			isItemPacked[i] += isItemPackedInCoR[i] + isItemPackedInCoN[i];
			model.addConstr(isItemPacked[i] == isItemPackedV[i]);
			// first quarter
			for (int j = 0; j < allo.items.size(); j++){
				if(i != j){
					model.addConstr(isItemPackedInCoX[i] - isItemPackedInCoX[j] + allo.items[i].wid <= allo.W * (1 - isItemPackedInCoN[i] + aij[i][j])); 
					model.addConstr(isItemPackedInCoX[j] - isItemPackedInCoX[i] + allo.items[j].wid <= allo.W * (1 - isItemPackedInCoN[j] + bij[i][j])); 
					model.addConstr(isItemPackedInCoY[i] - isItemPackedInCoY[j] + allo.items[i].len <= allo.L * (1 - isItemPackedInCoN[i] + cij[i][j])); 
					model.addConstr(isItemPackedInCoY[j] - isItemPackedInCoY[i] + allo.items[j].len <= allo.L * (1 - isItemPackedInCoN[j] + dij[i][j])); 
					model.addConstr(aij[i][j] + bij[i][j] + cij[i][j] + dij[i][j] <= 3); 
				}
			}
			// second quarter
			for (int j = 0; j < allo.items.size(); j++){
				model.addConstr(isItemPackedInCoX[i] - isItemPackedInCoX[j + allo.nbItems] + allo.items[i].wid <= allo.W * (1 - isItemPackedInCoN[i] + aij[i][j+allo.nbItems])); 
				model.addConstr(isItemPackedInCoX[j + allo.nbItems] - isItemPackedInCoX[i] + allo.items[j].widR <= allo.W * (1 - isItemPackedInCoR[j] + bij[i][j + allo.nbItems])); 
				model.addConstr(isItemPackedInCoY[i] - isItemPackedInCoY[j + allo.nbItems] + allo.items[i].len <= allo.L * (1 - isItemPackedInCoN[i] + cij[i][j+allo.nbItems])); 
				model.addConstr(isItemPackedInCoY[j+ allo.nbItems] - isItemPackedInCoY[i] + allo.items[j].lenR <= allo.L * (1 - isItemPackedInCoR[j] + dij[i][j+allo.nbItems])); 
				model.addConstr(aij[i][j + allo.nbItems] + bij[i][j + allo.nbItems] + cij[i][j + allo.nbItems] + dij[i][j + allo.nbItems] <= 3);
			}
			// third quarter
			for (int j = 0; j < allo.items.size(); j++){
				model.addConstr(isItemPackedInCoX[i+ allo.nbItems] - isItemPackedInCoX[j] + allo.items[i].widR <= allo.W * (1 - isItemPackedInCoR[i] + aij[i+ allo.nbItems][j])); 
				model.addConstr(isItemPackedInCoX[j] - isItemPackedInCoX[i+ allo.nbItems] + allo.items[j].wid <= allo.W * (1 - isItemPackedInCoN[j] + bij[i+ allo.nbItems][j])); 
				model.addConstr(isItemPackedInCoY[i+ allo.nbItems] - isItemPackedInCoY[j] + allo.items[i].lenR <= allo.L * (1 - isItemPackedInCoR[i] + cij[i+ allo.nbItems][j])); 
				model.addConstr(isItemPackedInCoY[j] - isItemPackedInCoY[i+ allo.nbItems] + allo.items[j].len <= allo.L * (1 - isItemPackedInCoN[j] + dij[i+ allo.nbItems][j])); 
				model.addConstr(aij[i+ allo.nbItems][j] + bij[i+ allo.nbItems][j] + cij[i+ allo.nbItems][j] + dij[i+ allo.nbItems][j] <= 3);
			}	
			// last quarter
			for (int j = 0; j < allo.items.size(); j++){
				if(i != j){
					model.addConstr(isItemPackedInCoX[i+ allo.nbItems] - isItemPackedInCoX[j+ allo.nbItems] + allo.items[i].widR <= allo.W * (1 - isItemPackedInCoR[i] + aij[i+ allo.nbItems][j+ allo.nbItems])); 
					model.addConstr(isItemPackedInCoX[j+ allo.nbItems] - isItemPackedInCoX[i+ allo.nbItems] + allo.items[j].widR <= allo.W * (1 - isItemPackedInCoR[j] + bij[i+ allo.nbItems][j+ allo.nbItems])); 
					model.addConstr(isItemPackedInCoY[i+ allo.nbItems] - isItemPackedInCoY[j+ allo.nbItems] + allo.items[i].lenR <= allo.L * (1 - isItemPackedInCoR[i] + cij[i+ allo.nbItems][j+ allo.nbItems])); 
					model.addConstr(isItemPackedInCoY[j+ allo.nbItems] - isItemPackedInCoY[i+ allo.nbItems] + allo.items[j].lenR <= allo.L * (1 - isItemPackedInCoR[j] + dij[i+ allo.nbItems][j+ allo.nbItems])); 
					model.addConstr(aij[i+ allo.nbItems][j + allo.nbItems] + bij[i+ allo.nbItems][j + allo.nbItems] + cij[i+ allo.nbItems][j + allo.nbItems] + dij[i+ allo.nbItems][j + allo.nbItems] <= 3);
				}
			}			
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
		
		// Global constraints -- linking expr and var
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
		
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit,  3600);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol,1e-9);
		model.optimize();
		
		// Filling Info
		allo.infos.timeCPU[0] = getCPUTime() - initTimeModelCPU;
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
		for (int i = 0; i < allo.items.size(); i++){	
			if(ceil(isItemPackedInCoN[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				cout << i << " " << isItemPackedInCoY[i].get(GRB_DoubleAttr_X) << " " << isItemPackedInCoX[i].get(GRB_DoubleAttr_X) << endl;
				allo.pps[i].push_back(ceil(isItemPackedInCoY[i].get(GRB_DoubleAttr_X)-EPSILON)); allo.pps[i].push_back(ceil(isItemPackedInCoX[i].get(GRB_DoubleAttr_X)-EPSILON)); allo.pps[i].push_back(0); 
				allo.columns[i].push_back(i); 
				for(int m = 1; m<allo.maxStack;m++){
					for (int l = 0; l < allo.items.size(); l++){
						if(allo.cAdj[i][l] && ceil(isItemStackedInCo[i][l][m].get(GRB_DoubleAttr_X) - EPSILON) == 1)
							allo.columns[i].push_back(l); 
					}
				}
			}
			if(ceil(isItemPackedInCoR[i].get(GRB_DoubleAttr_X) - EPSILON) == 1){
				cout << i << " " << isItemPackedInCoY[i + allo.nbItems].get(GRB_DoubleAttr_X) << " " << isItemPackedInCoX[i + allo.nbItems].get(GRB_DoubleAttr_X) << endl;
				allo.pps[i].push_back(ceil(isItemPackedInCoY[i+allo.nbItems].get(GRB_DoubleAttr_X)-EPSILON)); allo.pps[i].push_back(ceil(isItemPackedInCoX[i+allo.nbItems].get(GRB_DoubleAttr_X)-EPSILON)); allo.pps[i].push_back(1); 
				allo.columns[i].push_back(i); 
				for(int m = 1; m<allo.maxStack;m++){
					for (int l = 0; l < allo.items.size(); l++){
						if(allo.cAdj[i][l] && ceil(isItemStackedInCo[i][l][m].get(GRB_DoubleAttr_X) - EPSILON) == 1)
							allo.columns[i].push_back(l); 
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