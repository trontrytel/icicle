#include <cgicc/Cgicc.h>
#include <cgicc/HTTPHTMLHeader.h>
#include <cgicc/HTMLClasses.h>

#include <unistd.h>
#include <cstring>
#include <sstream>

#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;
using namespace cgicc;

#define error_macro(msg) \
{ \
  cout << cgicc::HTTPHTMLHeader() << msg << endl; \
  exit(EXIT_SUCCESS); \
}

void refresh(const string &msg)
{
  cout 
    << HTTPHTMLHeader() << endl
    << html()
    << head()
    << meta()
      .set("http-equiv","refresh")
      .set("content","10")
    << script()
      .set("src", "http://fgnass.github.com/spin.js/dist/spin.min.js")
      .set("type","text/javascript")
    << script()
    << head()
    << body()
      .set("style","text-align:center") 
    << cgicc::div()
      .set("id","spin")
      .set("style","padding:8em") 
    << cgicc::div()
    << script()
      .set("type","text/javascript")
    << "var target = document.getElementById('spin');" << endl
    << "var spinner = new Spinner().spin(target);" << endl
    << script()
    << h3() << msg << h3() << endl
    << pre().set("style","font-size:1.5em; margin-top: 5em;");
  cout.flush();
  int status = system("/usr/games/fortune -s");
  cout 
    << pre()
    << body()
    << html();
  exit(EXIT_SUCCESS);
}

int main(int argc, char **argv)
{
  Cgicc cgi;

  char *qs = getenv("QUERY_STRING");

  if (string(qs) == "purge") 
  {
    int status = system("rm -rf http/*");
    error_macro("purged.");
  }

  //if (qs == NULL) abort("QUERY_STRING undefined");
  boost::hash<string> hash;

  ostringstream tmp;
  tmp << "http/" << hash(qs);
  string dir = tmp.str();

  if (!boost::filesystem::exists(dir))
  {
    string cmd = 
      "touch " + dir + ".running;" + 
      "./test_todo_calc " + dir + ";" +
      "./test_todo_plot " + dir + ";" +
      "rm " + dir + ".running &";

    int pid = fork();
    if (pid < 0) error_macro("fork() failed.")
    else if (pid != 0) refresh("simulation started, please wait...");

    FILE* status;
    status = freopen( "/dev/null", "r", stdin); 
    status = freopen( "/dev/null", "w", stdout);
    //freopen( "/dev/null", "w", stderr); // keeping to log error messages
    execl("/bin/bash", "bash", "-c", cmd.c_str(), NULL);
  }
  else
  {
    if (!boost::filesystem::exists(dir + ".running")) 
      cout << "Location: " << dir << "/todo.gif" << endl << endl;
    else 
      refresh("please wait a bit longer...");
  }
}

