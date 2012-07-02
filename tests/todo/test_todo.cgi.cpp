#include <cgicc/Cgicc.h>
#include <cgicc/HTTPHTMLHeader.h>

#include <unistd.h>
#include <cstring>
#include <sstream>

#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>

#define error_macro(msg) { \
  std::cout << cgicc::HTTPHTMLHeader() << std::endl; \
  std::cout << msg << std::endl; \
}

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  cgicc::Cgicc cgi;

  char *qs = getenv("QUERY_STRING");
  if (qs == NULL) error_macro("QUERY_STRING undefined");
  boost::hash<string> hash;

  ostringstream tmp;
  tmp << "http/" << hash(qs);
  string dir = tmp.str();

  if (!boost::filesystem::exists(dir))
  {
    string cmd = "./test_todo_calc " + dir + " && ./test_todo_plot " + dir + " &";
    int status = system(cmd.c_str()); 
    if (status == -1) error_macro(cmd << " said: " << strerror(errno)) // TODO: errno
    error_macro("simulation running...")
  }
  else
  {
    if (boost::filesystem::exists(dir + "/test_00000.png")) 
      cout << "Location: " << dir << "/test_00001.png" << endl << endl;
    else error_macro("please refresh...");
  }
}

