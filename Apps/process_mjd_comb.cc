#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <bitset>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 4) {
    cout << "Usage: " << argv[0] << " [file gat ] [file enr] [outputfile]" << endl;
    return 1;
  }
  
  string file1 = argv[1];
  string file2 = argv[2]; 
  string file3 = argv[3];
  /*
  TFile *f1 = new TFile(Form("%s",file1.c_str()));
  TTree *t1 = (TTree*)f1->Get("mjdTree");
  TFile *f2 = new TFile(Form("%s",file2.c_str()));
  TTree *t2 = (TTree*)f2->Get("mjdTree");
  //TChain* f1 = new TChain("gatTree");
  //TChain* f2 = new TChain("enrTree");
  //f1->AddFile(file1.c_str());
  //f2->AddFile(file2.c_str());
  TFile *file = new TFile(file3.c_str(),"RECREATE");
  t1->AddFriend("mjdTree",f2);
  t2->AddFriend("mjdTree",f1);
  t1->Fill();
  t1->Print();
  //t2->Print();
  //f2->AddFriend(f1);
  TList *list1 = t1->GetListOfFriends();
  //list1->Print();
  */
  TChain ch("mjdTree");
  ch.Add(file1.c_str());
  ch.Add(file2.c_str());
  ch.Merge(file3.c_str());
  //  TFile *file = new TFile(file3.c_str(),"RECREATE");
  
  //TTree *newTree = t1->CloneTree(0);
  //TTree *newTree = t1->CopyTree();
  //TTree *newTree = t1->CloneTree(-1,"fast");
  //newTree->Print();
  //file->Write();
  //delete newTree;
  //delete file;
  //delete f1;
  //delete f2;
}


