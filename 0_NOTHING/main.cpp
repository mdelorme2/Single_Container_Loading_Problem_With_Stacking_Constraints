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
		vector<vector<vector<GRBVar> > > isItemPackedInCoN (allo.nbItems);
		vector<vector<vector<GRBVar> > > isItemPackedInCoR (allo.nbItems);
		vector<vector<vector<GRBVar> > > isItemStackedInCo (allo.nbItems);
		vector<vector<GRBLinExpr> > isSpaceBusy(allo.L);
		vector<GRBLinExpr> isItemStacked(allo.nbItems,0);
		vector<GRBLinExpr> isItemPacked(allo.nbItems,0);
		vector<GRBLinExpr> heighColumn(allo.nbItems,0);
		vector<GRBVar> isItemStackedV(allo.nbItems);
		vector<GRBVar> isItemPackedV(allo.nbItems);

		// Initialization
		for (int i = 0; i < allo.items.size(); i++){	
			isItemStackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedV[i] = model.addVar(0, 1, 0, GRB_BINARY);
			isItemPackedInCoN[i].resize(allo.L);
			isItemPackedInCoR[i].resize(allo.L);
			for (int j = 0; j < allo.L; j++){	
				if(allo.items[i].len <= allo.L){
					isItemPackedInCoN[i][j].resize(allo.W);
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].wid <= allo.W)
							isItemPackedInCoN[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);
					}
				}
				if(j + allo.items[i].lenR <= allo.L){
					isItemPackedInCoR[i][j].resize(allo.W);
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].widR <= allo.W)
							isItemPackedInCoR[i][j][k] =  model.addVar(0, 1, 0, GRB_BINARY);								
					}
				}
			}
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
		for (int i = 0; i < allo.items.size(); i++){	
			for (int j = 0; j < allo.L; j++){	
				if(j + allo.items[i].len <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].wid <= allo.W){
							isItemPacked[i] += isItemPackedInCoN[i][j][k];
							for(int j2 = 0; j2 < allo.items[i].len ; j2++){
								for(int k2 = 0; k2 < allo.items[i].wid ; k2++)
									isSpaceBusy[j+j2][k+k2] += isItemPackedInCoN[i][j][k];
							}
						}
					}
				}
				if(j + allo.items[i].lenR <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].widR <= allo.W){
							isItemPacked[i] += isItemPackedInCoR[i][j][k];
							for(int j2 = 0; j2 < allo.items[i].lenR ; j2++){
								for(int k2 = 0; k2 < allo.items[i].widR ; k2++)
									isSpaceBusy[j+j2][k+k2] += isItemPackedInCoR[i][j][k];
							}
						}
					}
				}
			}
			model.addConstr(isItemPacked[i] == isItemPackedV[i]);
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
		
		// Nooverlap
		for (int k = 0; k < allo.L; k++){	
			for (int l = 0; l < allo.W; l++)
				model.addConstr(isSpaceBusy[k][l] <= 1);
		}
		
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
			for (int j = 0; j < allo.L; j++){	
				if(j + allo.items[i].len <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].wid <= allo.W && ceil(isItemPackedInCoN[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
							allo.pps[i].push_back(j); allo.pps[i].push_back(k); allo.pps[i].push_back(0); 
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
			}
			for (int j = 0; j < allo.L; j++){	
				if(j + allo.items[i].lenR <= allo.L){
					for (int k = 0; k < allo.W; k++){
						if(k + allo.items[i].widR <= allo.W && ceil(isItemPackedInCoR[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
							allo.pps[i].push_back(j); allo.pps[i].push_back(k); allo.pps[i].push_back(1); 
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