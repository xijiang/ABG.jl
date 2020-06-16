#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <thread>
#include <utility>
#include <cmath>
#include <numeric>
#include <vector>
#include <tuple>
#include <filesystem>
#include <openblas/cblas.h>

/**
 * This program is designed only for Linux system.
 * To use multiple threads, it is better recompile openblas multiple-threaded
 * Then,
 * set OPENBLAS_NUM_THREADS environment variable to specify threads to be used.
 *
 * To use the current version:
 * g++ -O2 -Wall -std=c++17 -lopenblas -lpthread bigg.cpp -o bigg
 */
using namespace std;

tuple<size_t, size_t, size_t> dims(istream&ins){
    // Determine number of loci, ID
    ins.seekg(0, ins.end);
    size_t fsz = ins.tellg();
    ins.seekg(0, ins.beg);
    
    string line;
    getline(ins, line);
    ins.seekg(0, ins.beg);
    size_t lsz = line.length() + 1; // length of each line, including '\n'
    size_t nlc = lsz - 1;
    size_t nid = fsz/lsz;
    
    clog << setw(36) << "Number of loci: " << setw(8) << nlc << endl;
    clog << setw(36) << "Number of ID: "  << setw(8) << nid << endl;
    if(lsz*nid != fsz) throw runtime_error("Error: Not a sqaure matrix.");
    return {nid, nlc, lsz};
}

size_t mem_plan(const size_t nid, const size_t nlc, double&mem, unsigned&nthreads){
  clog << fixed;
  clog.precision(2);
  double gib = 1024 * 1024 * 1024;

  // CPU related
  unsigned tth = thread::hardware_concurrency();
  clog << setw(36) << "Number of CPU threads: " << setw(8) << tth << endl;
  clog << setw(36) << "Number of threads by default: " << setw(8) << nthreads << endl;
  if(nthreads > tth) nthreads = tth-1;
  if(!nthreads) nthreads = 1;
  clog << setw(36) << "Number of threads assigned: " << setw(8) << nthreads << endl;
  
  // Memory related
  ifstream fin("/proc/meminfo");{
    string line;
    getline(fin, line);
    stringstream ts(line);
    double tm;
    ts>>line>>tm;
    tm = tm*1024/gib;
    clog << setw(36) << "Total system memory: " << setw(8) << tm << " GiB" << endl;
    getline(fin, line);
    stringstream fs(line);
    fs>>line>>tm;
    tm = tm*1024/gib;
    clog << setw(36) << "Free memory: " << setw(8) << tm << " GiB" << endl;
    clog << setw(36) << "Memory asked: " << setw(8) << mem << " GiB" << endl;
    if(mem > tm) mem = tm;
    clog << setw(36) << "Memory allowed: " << setw(8) << mem << " GiB" << endl;
  }
  
  // Block size
  double c = mem * gib/ 8;
  double sol = sqrt(nlc * nlc + c) - nlc;
  size_t nln = ceil(sol);
  if(nln >= nid) nln = nid;
  clog << setw(36) << "Number of lines a time: " << setw(8) << nln << endl;
  
  // Disk spaces
  double gmt = 4. * nid * (nid-1) / gib;
  clog << setw(36) << "Disk required for final G matrix: ";
  clog << setw(8) << gmt << " GiB" << endl;
  double mid{0.};
  mid  = gmt;
  mid += nid/nln * (nln-1)*(nln-2)/2 * 8 / gib;
  mid += (nid % nln - 1) * (nid % nln - 2) / 2 * 8 /gib;

  clog << setw(36) << "Disk usage for mid-results: "
       << setw(8) << mid << " GiB" << endl;
  filesystem::space_info pwd = filesystem::space(".");
  double avail = pwd.available/gib;
  clog << setw(36) << "Disk space available: " << setw(8) << avail << " GiB" << endl;
  if(avail - .1 < gmt + mid) throw runtime_error("Not enough disk space");
  return nln;
}

void read_mblk(istream&ins, vector<double>&m, const size_t&lsz,
	       size_t beg, size_t end, const vector<double>&twop){
  size_t k{0};
  string line;
  ins.seekg(lsz * beg);
  for(size_t i=beg; i<end; ++i){
    getline(ins, line);		// prepare for vanRaden method I
    for(size_t j=0; j<line.length(); ++j) m[k++] = line[j] - '0' - twop[j];
  }
}

void title(const string&msg){
  clog<<endl<<msg<<endl;
  for(size_t i=0; i<msg.length()+2; ++i) clog<<'=';
  clog << endl;
}

void nuclear(){
  return;
}

void calc_block(istream&raw,
		const unsigned&nth,
		const vector<double>&frq,
		const size_t&lsz,
		const vector<pair<size_t, size_t>>&range,
		const vector<pair<int, int>>&tasks){ // vanRaden method I
  //double alpha{1.};
  //double beta{0.};
  thread tsk[nth];
  //size_t i, j, k, nblk(range.size()), nlc(lsz-1);
  unsigned t{0};
  for(auto&[x,y]:tasks){
    tsk[t++] = thread(nuclear);
    if(t == nth) for(unsigned i{0}; i<nth; ++i) tsk[i].join();
  }
  if(t) for(unsigned i{0}; i<t; ++i) tsk[i].join();
  /*
  // the nuclear computation block
  for(i=0; i<nblk; ++i){
    clog << '\r' << setw(15) << "block row:" << setw(3) << i+1
	 << '/' << nblk << ';' << "\tcolumn: " << setw(3) << 1
	 << '/' << nblk-i << "  " << flush;
    read_mblk(raw, mi, lsz, range[i].first, range[i].second, twop);
    size_t m = range[i].second - range[i].first;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, nlc,
		alpha, &mi[0], nlc, &mi[0], nlc, beta, &gblk[0], m);
    for(k=0; k<m*m; ++k) gblk[k]/=s2pq;
    ofstream fii(to_string(i) + ".g");
    fii.write(reinterpret_cast<char*>(&gblk[0]), sizeof(double)*m*m);
    for(j = i+1; j<nblk; ++j){
      clog << '\r' << setw(15) << "block row:" << setw(3) << i+1
	   << '/' << nblk << ';' << "\tcolumn: " << setw(3) << j-i+1
	   << '/' << nblk-i << "  " << flush;

      size_t n = range[j].second - range[j].first;
      read_mblk(raw, mj, lsz, range[j].first, range[j].second, twop);
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, nlc,
		  alpha, &mi[0], nlc, &mj[0], nlc, beta, &gblk[0], n);
      for(k=0; k<m*n; ++k) gblk[k]/=s2pq;
      ofstream fij(to_string(i)+'-'+to_string(j)+".g");
      fij.write(reinterpret_cast<char*>(&gblk[0]), sizeof(double)*m*n);
    }
  }
  clog << endl;
  */
}

void merge_block(const vector<pair<size_t, size_t>>&range){
  size_t i, j, m, n, nblk(range.size());
  for(i=0; i<nblk; ++i){
    m = range[i].second - range[i].first;
    vector<double> diag(m*m);
    ifstream fii(to_string(i) + ".g");
    fii.read(reinterpret_cast<char*>(&diag[0]), sizeof(double)*m*m);
    for(j=i+1; j<nblk; ++j){
      n = range[i].second - range[j].first;
      vector<double> offd(m*n);
      ifstream fij(to_string(i)+'-'+to_string(j)+".g");
      fij.read(reinterpret_cast<char*>(&offd[0]), sizeof(double)*m*n);
    }
  }
}

int main(int argc, char *argv[])
{
  ios_base::sync_with_stdio(false);
  
  if(argc != 3){
    cerr << "Usage: cat freq-file | " << argv[0] << " genotype-file memory\n";
    cerr << "memory are the amount of required RAM in GiB\n";
    cerr << "this program uses 6 threads by default\n";
    return 1;
  }
  
  vector<double> frq;		// Read frequencies from stdin.
  for(double f; cin>>f; frq.push_back(f));
  
  ifstream raw(argv[1]);
  if(raw){			// file successfully opened.
    ////////////////////////////////////////////////////////////
    title("System summary");
    // Determine parameters.
    auto [nid, nlc, lsz] = dims(raw);
    double mem(stof(argv[2]));
    unsigned nthreads(6);
    auto nln = mem_plan(nid, nlc, mem, nthreads);
    if(nlc != frq.size()){
      throw runtime_error("Numbers of loci don't agree");
      return 2;
    }
    //vector<double> gblk(nln * nln), mi(nln * nlc), mj(nln * nlc);

    ////////////////////////////////////////////////////////////
    title("The calculation procedure");
    // debug--begin
    nln = 1;
    // debug--end
    size_t nblk(nid/nln);
    if(nid%nln) ++nblk;
    clog << setw(36) << "The genotypes will be dealt with in "<<nblk<<" blocks"<<endl;
    vector<pair<size_t, size_t>> range;
    vector<pair<int, int>> tasks;
    for(size_t i=0; i<nblk; ++i) range.push_back({i*nln, i*nln+nln});
    //range[nblk-1].second = nid; // blocks determined.
    //for(auto i{0}; i<static_cast<int>(range.size()); ++i)
    //  for(auto j{i}; j<static_cast<int>(range.size()); ++j)
    //	tasks.push_back({i, j});
    //calc_block(raw, nthreads, frq, lsz, range, tasks);
    //
    //title("Merge blocks");
    //merge_block(range);
    // them merge the small blocks into one big-G
  }else{
    cerr << "Invalid file: " << argv[1] <<'\n';
    return 3;
  }

  return 0;
}
