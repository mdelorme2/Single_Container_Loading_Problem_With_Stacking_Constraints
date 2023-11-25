#include "Allocation.h" 

/*	*************************************************************************************
	************************************* ITEM ******************************************
	************************************************************************************* */

void Item::print(){
	cout << idx << " " << id << " " << col << " " << lenO << " " << widO << " " << heiO << " " << pro << " " << wei << " " << sta << " " << fla  << " [" << len << " " << wid  << " " << hei << "] [" << lenR << " " << widR << " " << hei << "]" << endl;
}
	
/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */
	

void Allocation::load(const string& path, const string& filein){
	// Local variables 
	istringstream iss;
	istringstream tempIss;
	string parser;
	int garbage;
	string nameFile = path + filein;
	string tempString;

	// File opening
	ifstream file(nameFile.c_str(), ios::in);

	// File lecture
	if (file){
		// Name of the instance is filein
		name = filein;

		// Truck info
		sumVol = 0;	
		getline(file, parser); iss.str(parser); 
		iss >> L; LO = L;
		iss >> W; WO = W;
		iss >> H; HO = H;
		iss >> maxStack;
		iss >> maxWeight;
		iss >> nbItems;
		iss.clear(); 
		
		for(int i = 0; i < nbItems; i++){
			Item it;
			getline(file, parser); iss.str(parser); 
			iss >> it.id; it.idx = i;
			iss >> it.len;
			iss >> it.wid;
			if(it.len < it.wid){
				int temp = it.len; it.len = it.wid; it.wid = temp;
			}
			iss >> it.hei;
			iss >> it.wei;
			iss >> it.sta;
			iss >> it.fla;
			iss >> it.col;
			it.lenO = it.len;
			it.widO = it.wid;
			it.lenR = it.wid;
			it.widR = it.len;
			it.heiO = it.hei;
			it.pro = it.lenO * it.widO /** it.heiO*/;
			sumVol += it.pro;
			items.push_back(it);
			iss.clear(); 
		}	
				
		file.close();	
	}
	else cout << "Could not open the file " << nameFile << endl;

	// Create compatibility graph
	cAdj.resize(items.size());
	for(int i =0; i<items.size();i++){
		cAdj[i].resize(items.size(),0);
		for(int j =0; j<items.size();j++){
			if (i != j && items[i].hei + items[j].hei <= H && items[i].sta == 1 && items[i].fla == 1 && items[j].sta == 1 && items[i].wei >= items[j].wei 
				&& (i < j || cAdj[j][i] == 0) && ((items[i].wid >= items[j].wid && items[i].len >= items[j].len) || (items[i].wid >= items[j].widR && items[i].len >= items[j].lenR))){
				if(items[i].col == items[j].col){
					cArcs.push_back({i,j});
					cAdj[i][j] = 1;
				}
			}
		}
	}	
}

void Allocation::printProb(){
	cout << "Instance " << name << endl;
	cout << "Truck dimension " << L << " x " << W << " x " << H << " - " << maxStack << " stacks - " <<  maxWeight << " kg in total "<< endl;
	for(int i = 0; i < items.size();i++){
		items[i].print(); 
	}
	cout << "Stackability" << endl;
	for(int i = 0; i < items.size();i++){
		cout << i << " -- ";
		for(int j = 0; j < items.size();j++){		
			if(cAdj[i][j])
				cout << j << " ";
		}
		cout << endl;
	}
	cout << "-----------------------------------------" << endl;
}

void Allocation::printInfo(const string& pathAndFileout){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << name << "\t" << infos.opt << "\t";
	for (int i = 0; i < infos.timeCPU.size();i++)
		file << infos.timeCPU[i] << "\t";
	file << infos.LB << "\t" << infos.UB <<  "\t"  << sumVol <<  "\t" << infos.nbVar << "\t" << infos.nbCons << "\t" << infos.nbNZ << "\t" << infos.nbCuts1 << "\t" << infos.nbCuts2 << "\t" << checkSolution() << endl;
	file.close();
}

void Allocation::printFile(const string& picture){
	string nameFile = picture;
	ofstream file(nameFile.c_str(), std::ios::out | ios::trunc);
	int nbColumns = 0;
	for(int i =0;i<items.size();i++){ 
		if(pps[i].size() > 0) 
			nbColumns++;
	}
	
	file << LO << " " << WO << " " << nbColumns << endl;
	for (int i = 0; i < items.size();i++){
		if(pps[i].size() > 0){
			file << pps[i][0] << " " << pps[i][1] << " " << items[i].col << " ";
			if(pps[i][2] == 0) file << items[i].lenO << " " << items[i].widO << " ";
			else file << items[i].widO << " " << items[i].lenO << " ";			
			for (int j = 0; j < columns[i].size();j++){
				// file << items[columns[i][j]].id << " ";
				file << columns[i][j] << " ";
			}
			file << endl;
		}
	}
	file.close();
}

void Allocation::rescaleW(){
	vector<int> trueW (W+1,0);
	bool change = true;
	while(change){
		change = false;
		for(int i = 0; i < items.size();i++){
			if(pps[i].size() == 0) continue;
			if(pps[i][2] == 0 && trueW[pps[i][1] + items[i].wid] < trueW[pps[i][1]] + items[i].widO){
				trueW[pps[i][1] + items[i].wid] = trueW[pps[i][1]] + items[i].widO;
				change = true;
			}
			if(pps[i][2] == 1 && trueW[pps[i][1] + items[i].widR] < trueW[pps[i][1]] + items[i].lenO){
				trueW[pps[i][1] + items[i].widR] = trueW[pps[i][1]] + items[i].lenO;
				change = true;
			}
			cout << "After dealing with " ;
			cout << i << " ";
			for (int j = 0; j < pps[i].size();j++){
				cout << pps[i][j] << " ";
			}		
			cout << "--- ";
			for(int j = 0; j < W + 1; j++){
				cout << trueW[j] << " ";
			}
			cout << endl;
		}	
		for(int i = 0; i < W ; i++){
			if(trueW[i+1] < trueW[i]){
				trueW[i+1] = trueW[i]; 
				change = true;
			}
		}
		cout << "After iteration ";
		for(int i = 0; i < W + 1; i++){
			cout << trueW[i] << " ";
		}
		cout << endl;
	}
	for(int i = 0; i < items.size();i++){
		if(pps[i].size() > 0)
			pps[i][1] = trueW[pps[i][1]];
	}
}
		
int Allocation::checkSolution(){
	bool isOkay = 1;
	vector<int> isItem(nbItems,0); 
	vector<vector<int> > isBusy(LO,vector<int> (WO,0));
	int weight = 0;
	int profit = 0;
	if(pps.size() == 0) return 0;
	
	// Start check
	for(int i =0;i<items.size();i++){ 
		if(pps[i].size() > 0){
			// Step 1 check that the column fits into the truck
			if(pps[i][2] == 0){	
				if(pps[i][0] + items[i].lenO > LO || pps[i][1] + items[i].widO > WO){
					cout << "Item "<< i << " does not fit in the truck " << endl;
					isOkay = 0;
				}
			}
			else{
				if(pps[i][0] + items[i].widO > LO || pps[i][1] + items[i].lenO > WO){
					cout << "Item "<< i << " does not fit in the truck " << endl;
					isOkay = 0;
				}
			}
			// Step 2 check that the column does not overlap with another column
			if(pps[i][2] == 0){	
				for(int j = pps[i][0]; j < pps[i][0] + items[i].lenO; j++){
					for(int k = pps[i][1]; k < pps[i][1] + items[i].widO; k++){
						if(isBusy[j][k] == 0) isBusy[j][k] = 1;
						else{
							cout << "Item "<< i << " overlaps with another item " << endl;
							isOkay = 0;
						}
					}
				}
			}
			else{
				for(int j = pps[i][0]; j < pps[i][0] + items[i].widO; j++){
					for(int k = pps[i][1]; k < pps[i][1] + items[i].lenO; k++){
						if(isBusy[j][k] == 0) isBusy[j][k] = 1;
						else{
							cout << "Item "<< i << " overlaps with another item " << endl;
							isOkay = 0;
						}
					}
				}
			}
			// Step 3 check that the column fits into the truck + that the column is feasible + compute total weight + compute total profit
			int height = 0;
			for (int j = 0; j < columns[i].size();j++){
				height += items[columns[i][j]].hei;
				profit += items[columns[i][j]].pro;
				weight += items[columns[i][j]].wei;
				if(isItem[columns[i][j]] == 0) isItem[columns[i][j]] = 1;
				else{
					cout << "Item "<< columns[i][j] << " is used more than once " << endl;
					isOkay = 0;
				}
				if(j >= 1 && cAdj[columns[i][j-1]][columns[i][j]] != 1){
					cout << "Item "<< columns[i][j] << " does not fit on top of item " << columns[i][j-1] << endl;
					isOkay = 0;					
				}
			}
			if(height > HO){
				cout << "Column "<< i << " is too high " << endl;
				isOkay = 0;
			}
			if( columns[i].size() > maxStack){
				cout << "Column "<< i << " has too many stacked items " << endl;
				isOkay = 0;					
			}
		}
	}
	
	// Step 4 check that the weight and profit are okay
	if(weight > maxWeight){
		cout << "Total weight above the threshold" << endl;
		isOkay = 0;			
	}
	if(profit != infos.LB){
		cout << "Total profit of the solution is not what it is supposed to be" << endl;
		isOkay = 0;				
	}
	return isOkay;
}