#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

void splitarray(vector<int>* inputarray, vector<int>* outputarray1,vector<int>* outputarray2){
  int size = inputarray->size();
  if(size%2 == 0){
    for(size_t i=0;i<inputarray->size();i++){
      if(i<inputarray->size()/2){
	outputarray1->push_back(inputarray->at(i));
      }else{
	outputarray2->push_back(inputarray->at(i));
      }
    }
  }else{
    for(size_t i=0;i<inputarray->size();i++){
      if(i<(inputarray->size()-1)/2){
	outputarray1->push_back(inputarray->at(i));
      }else if(i>(inputarray->size()-1)/2){
	outputarray2->push_back(inputarray->at(i));
      } 
    }
  }
}

double median(vector<int>* inputarray){
  double mean = 0;
  int size = inputarray->size();  
  if(size%2==1){
    mean = inputarray->at((size-1)/2);
  }else{
    mean = (inputarray->at(size/2-1)+inputarray->at(size/2))/2;
  }
  return mean;
}

int main() {
  /* Enter your code here. Read input from STDIN. Print output to STDOUT */   
  int N;
  vector<int> array;
  //cout << "Please Input Array Size:";
  cin >> N;
  //cout << N << endl;
  //cout << "Please Input Array Item:";
  for(int i =0;i<N;i++){
    int j;
    cin >> j;
    //cout << "Item " << i << ": " << j << endl ;
    array.push_back(j);
  }
  
  //int size = array.size();

  std::sort(array.begin(),array.end());
  vector<int> output1;
  vector<int> output2;

  splitarray(&array, &output1,&output2);

  double quartile1 = median(&output1);
  double quartile2 = median(&array);
  double quartile3 = median(&output2);
  cout << quartile1 << endl;
  cout << quartile2 << endl;
  cout << quartile3 << endl;

  return 0;
}
