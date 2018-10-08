#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <ctype.h>
#include <vector>
#include <stdio.h> 
#include <algorithm>

using namespace std;


///////////////////////////////////////////////////////////////////////////////
//////////////////////////     IMPLEMENTATION     /////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void addFunc1(vector<pair<string,int>* >* toAdd, vector<pair<string,int> > toSort){

	for(int i=0; i<toAdd->size(); i++){
		pair<string,int>* currFinal = toAdd->at(i);
		for(int k=0; k<toSort.size(); k++){
			pair<string,int> currTemp = toSort.at(k);
			if(currFinal->first == currTemp.first){
				currFinal->second += 1;
			}
		}
	}
	
}

void addFunc(vector<pair<string,int>* >* toAdd){
	pair<string,int>* pluto = new pair<string, int>;

	pluto->first = "pluto";
	pluto->second = 9; 

	if(find(toAdd->begin(), toAdd->end(), pluto) == toAdd->end()) {
		toAdd->push_back(pluto);
	}
}

void modifyVec(vector<pair<string,int>* >* toModify){
	for(int i=0; i<toModify->size(); i++){
		pair<string,int>* curr = toModify->at(i);
		curr->second = curr->second + 1;
	}
}

vector<pair<string,int>* >* makeVec(){
	vector<pair<string,int>* >* toReturn = new vector<pair<string,int>* >; 

	pair<string,int>* mercury = new pair<string, int>;
	pair<string,int>* venus = new pair<string, int>;
	pair<string,int>* earth = new pair<string, int>;
	pair<string,int>* mars = new pair<string, int>;
	pair<string,int>* jupiter = new pair<string, int>;
	pair<string,int>* saturn = new pair<string, int>;
	pair<string,int>* uranus = new pair<string, int>;
	pair<string,int>* neptune = new pair<string, int>;

	mercury->first = "mercury";
	mercury->second = 1; 
	venus->first = "venus";
	venus->second = 2; 
	earth->first = "earth";
	earth->second = 3; 
	mars->first = "mars";
	mars->second = 4; 
	jupiter->first = "jupiter";
	jupiter->second = 5; 
	saturn->first = "saturn";
	saturn->second = 6; 
	neptune->first = "neptune";
	neptune->second = 7; 
	uranus->first = "uranus";
	uranus->second = 8; 
	
	toReturn->push_back(mercury); 
	toReturn->push_back(venus); 
	toReturn->push_back(earth); 
	toReturn->push_back(mars); 
	toReturn->push_back(jupiter); 
	toReturn->push_back(saturn); 
	toReturn->push_back(neptune); 
	toReturn->push_back(uranus); 


	return toReturn; 
}

int main() {
	cout << endl; 
	vector<pair<string,int>* >* myVec = makeVec();
	for(int i=0; i<myVec->size(); i++){
		pair<string,int>* curr = myVec->at(i);
		cout << curr->first << "   " << curr->second << endl;
	}
	cout << "------------------" << endl; 

	modifyVec(myVec); 
	for(int k=0; k<myVec->size(); k++){
		pair<string,int>* curr = myVec->at(k);
		cout << curr->first << "   " << curr->second << endl;
	}
	cout << "------------------" << endl; 

	addFunc(myVec); 
	for(int h=0; h<myVec->size(); h++){
		pair<string,int>* curr = myVec->at(h);
		cout << curr->first << "   " << curr->second << endl;
	}
	cout << "------------------" << endl << endl; 


	vector<pair<string,int> > smallList; 
	pair<string,int> mars, venus; 
	mars = make_pair("mars", 5);
	venus = make_pair("venus", 6);
	smallList.push_back(mars);
	smallList.push_back(venus);

	addFunc1(myVec, smallList);
	for(int l=0; l<myVec->size(); l++){
		pair<string,int>* curr = myVec->at(l);
		cout << curr->first << "   " << curr->second << endl;
	}
	cout << "------------------" << endl << endl; 
}



