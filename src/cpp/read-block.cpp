#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  if(argc != 5){
    cerr << "Usage:\n\t" << argv[0]  << " gt-file line-from line-to line-size\n\n";
    cerr << "\t\t- Line number start from 1\n";
    cerr << "\t\t- Line-from and -to inclusive\n";
    cerr << "\t\t- Line size include the newline character\n";
    cerr << "\t\t- Suppose newline is only one character long\n";
    return 1;
  }
  ifstream fin(argv[1]);
  size_t fra(stoi(argv[2]) - 1); // the from line number
  size_t til(stoi(argv[3]));	 // the to line number, exclusive
  size_t lsz(stoi(argv[4]));
  size_t nlc(lsz - 1), i{0}, k;
  size_t nid(til - fra);
  vector<double> dat(nid*nlc);
  string line;
  
  fin.seekg(fra*lsz);
  for(k = fra; k<til; ++k){
    getline(fin, line);
    for(auto&x:line) dat[i++] = x - '0';
  }
  cout.write(reinterpret_cast<char*>(&dat[0]), sizeof(double) * dat.size());
  return 0;
}
