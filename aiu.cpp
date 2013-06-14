#include <iostream>
#include <sstream>

using namespace std;


int main(void)
{
  ostringstream oss;
  oss << "aiu" << 10 << "aiu" << (long)100;
  cout << oss.good() << oss.eof() << oss.fail() << oss.bad() << endl;
  cout << oss.str() << endl;
  return 0;
}
