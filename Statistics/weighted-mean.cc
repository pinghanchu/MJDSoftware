#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;


int main() {
  /* Enter your code here. Read input from STDIN. Print output to STDOUT */   
  int N;
  vector<int> X;
  vector<int> W;
  //cout << "Please Input Array Size:";
  cin >> N;
  //cout << N << endl;

  for(int i =0;i<N;i++){
    int j;
    cin >> j;
    cout << j << endl;
    X.push_back(j);
  }
  for(int i =0;i<N;i++){
    int j;
    cin >> j;
    cout << j << endl;
    W.push_back(j);
  }
  double sum = 0;
  double count = 0;
  for(size_t ia = 0;ia<X.size();ia++){
    sum = (double)(sum+X.at(ia)*W.at(ia));
    count = (double)(count+W.at(ia));
  }

  double mean = sum/count;
  cout << setiosflags(ios::fixed) << setprecision(1) << mean << endl;
  return 0;
}
