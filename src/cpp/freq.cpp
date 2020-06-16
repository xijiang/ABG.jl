#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
  string line;
  getline(cin, line);
  int nlc = line.length(), nid{1};
  double frq[nlc];
  for(auto i{0}; i<nlc; ++i) frq[i] = line[i] - '0';
  while(getline(cin, line)){
    for(auto i{0}; i<nlc; ++i) frq[i] += line[i] - '0';
    ++nid;
  }
  nid *= 2;			// number of alleles now
  for(auto i{0}; i<nlc; ++i) cout << frq[i]/nid << '\n';
  return 0;
}
