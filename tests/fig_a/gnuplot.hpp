#pragma once

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

template <class data_t>
void plot(Gnuplot &gp, data_t data, std::string title)
{
  blitz::Array<float, 2> tmp(data);

  gp << "set title " << title << "\n";  
  gp << "splot '-' binary" << gp.binfmt(tmp) /*<< dxdy*/ << " with image failsafe notitle\n";

  gp.sendBinary(tmp);
}
