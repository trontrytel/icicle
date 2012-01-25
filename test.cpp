#include <blitz/array.h>
using namespace blitz;
int main()
{
  Array<float, 1> a(Range(0,10));
  a(11) = 1;
}
