/**
 * This program read pedigree from stdin.
 * It then calculate inbreeding coefficient of everybody.
 */

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <limits>

using namespace std;
using PM =pair<int, int>;
using PED=vector<PM>;
using MID=map<PM, double>;	// store intermediate resultes

double Amat(int i, int j, const PED &ped, MID&mid){
  if(i==0 || j==0) return 0;	// so that relationship with un unknown is not stored in mid
  
  if(i>j) swap(i, j);		// in <algorithms>, but may be included by utility or limits.
  
  // Look up if {i,j} was calculated before
  if(mid.find({i, j}) != mid.end()) return mid[{i, j}];

  const auto &[pa, ma] = ped[j];
  if(i==j)
    mid[{j, j}] = 1 + Amat(pa, ma, ped, mid) / 2.;
  else
    mid[{i, j}] = (Amat(i, pa, ped, mid) + Amat(i, ma, ped, mid)) / 2.;

  return mid[{i, j}];
}

int main(int argc, char *argv[])
{
  PED ped;			// the pedigree to look up
  int nid{0};

  ped.push_back({0,0});	// the magic dummy
  // hence, row number is ID number id starts from 1.
  for(int pa, ma; cin>>pa>>ma; ped.push_back({pa, ma})){
    ++nid;
    if(pa >= nid ||
       ma >= nid ||
       pa <  0   ||
       ma <  0      ){		// error check on the pedigree
      cerr << "ERROR: Invalid pa / ma ID @ line number: " << nid << endl;
      return 1;
    }
  }

  map<PM, double> mid;		// store the mid results of Amat
  
  typedef numeric_limits< double > dbl; // to avoid binary I/O
  cout.precision(dbl::max_digits10);

  for(auto id{1}; id<=nid; ++id)
    cout << Amat(id, id, ped, mid) - 1. <<'\n';

  return 0;
}
