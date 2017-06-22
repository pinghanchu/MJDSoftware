#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;


int main() {
  /* Enter your code here. Read input from STDIN. Print output to STDOUT */   
  int N;
  vector<int> array;
  //cout << "Please Input Array Size:";
  cin >> N;
  cout << N << endl;
  //cout << "Please Input Array Item:";
  for(int i =0;i<N;i++){
    int j;
    cin >> j;
    //cout << "Item " << i << ": " << j << endl ;
    array.push_back(j);
  }
  double mean = 0;
  double median = 0;
  double min = 1e6;
  double sum = 0;
  double count=0;
  int size = array.size();

  std::sort(array.begin(),array.end());
  for(size_t ia = 0;ia<array.size();ia++){
    sum = sum + (double)array.at(ia);
    count = count+1;
    if(array.at(ia)<min){
      min = array.at(ia);
    }
  }
  vector<int> item;
  vector<int> itemmod;
  item.push_back(array.at(0));
  itemmod.push_back(1);
  
  for(size_t ia=1;ia<array.size();ia++){
    if(array.at(ia)!=array.at(ia-1)){
      item.push_back(array.at(ia));
      itemmod.push_back(1);
    }else{
      itemmod.at(itemmod.size()-1)= itemmod.at(itemmod.size()-1)+1;
    }        
  }
  int itemmodmax = 0;
  for(size_t ia =0;ia<itemmod.size();ia++){
    if(itemmod.at(ia)>itemmodmax){
      itemmodmax = itemmod.at(ia);
    }
  }
  vector<int> itemmax;
  for(size_t ia = 0;ia<item.size();ia++){
    if(itemmod.at(ia)==itemmodmax){
      itemmax.push_back(item.at(ia));
    }
  }
  min = 1e6;
  for(size_t ia=0;ia<itemmax.size();ia++){
    if(itemmax.at(ia)<min){
      min = itemmax.at(ia);
    }
  }
  mean = sum/count;
  if(array.size()%2==0){
    median = (array.at(size/2-1)+array.at(size/2))/2.;
  }else{
    median = array.at(size/2);
  }
  cout << mean << endl;
  cout << median << endl;
  cout << min << endl;
  return 0;
}
