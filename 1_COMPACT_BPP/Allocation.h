#ifndef ALLOCATION_H
	#define ALLOCATION_H
	
	using namespace std;
	#include <iostream> 
	#include <iomanip> 
	#include <fstream>
	#include <sstream>
	#include <vector>
	#include <string>
	#include <set>
	#include <iostream> 
	#include <math.h> 
	#include <cstdlib>
	#include <algorithm>
	
	class Info;
	class Item;
	class Allocation;

/*	*************************************************************************************
	************************************* ITEM ******************************************
	************************************************************************************* */

	class Item{
	public:
		int id;
		int idx;
		int len;
		int lenR;
		int lenO;
		int wid;
		int widR;
		int widO;
		int hei;
		int heiO;
		int pro;
		int wei;
		int sta;
		int fla;
		int col;
		void print();
	};
	
/*	*************************************************************************************
	************************************* INFO ******************************************
	************************************************************************************* */

	class Info{
	public:
		bool opt;
		vector<double> timeCPU;
		int LB;
		int UB;
		float contUB;
		int nbCons;
		int nbVar;
		int nbNZ;
	};

/*	*************************************************************************************
	********************************** ALLOCATION ***************************************
	************************************************************************************* */

	class Allocation{
	public:
		
		// Data read from the file
		string name;
		int L; int LO;
		int W; int WO;
		int H; int HO;
		int maxStack;
		int maxWeight;
		int nbItems;
		vector<Item> items;
		int sumVol;
		
		// Necessary for the models
		vector<vector<int> > cArcs;
		vector<vector<int> > cAdj;
	
		// Given by the ILP model
		Info infos;
		vector<vector<int> > columns;
		vector<vector<int> > pps; 

		void load(const string& path, const string& filein);
		void printProb();
		void printInfo(const string& pathAndFileout);
		void printFile(const string& picture);
		void rescaleW();
		int checkSolution();
	};
	
	
#endif 