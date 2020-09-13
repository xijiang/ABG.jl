#include <iostream>
#include <filesystem>
/**
 * take the path from argv[1], then print the free and available space in bytes of current path.
 */

int main(int argc, char *argv[])
{
  if(argc!=2){
    std::cerr << "Error: path not specified\n";
    std::cout<<"0 0\n";
    return 1;
  }
  
  std::filesystem::space_info path = std::filesystem::space(argv[1]);
  std::cout <<path.free << ' ' << path.available << '\n';
  return 0;
}
