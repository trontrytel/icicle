#pragma once

#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

void init(Gnuplot &gp, const std::string &file, const int &ny, const int &nx)
{
  gp << "set term svg dynamic enhanced fsize 13 size " << nx * 500 << "," << ny * 500 << "\n";
  gp << "set size square\n";
  gp << "set encoding utf8\n";
  // progressive-rock connoisseur palette ;)
  gp << "set palette defined (0 '#FFFFFF', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n";
  gp << "set view map\n";
  gp << "dx = dy = 1500./75\n"; // TODO: read it from the HDF file
  gp << "set xtics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dx, '0.6' 600/dx, '0.9' 900/dx, '1.2' 1200/dx, '1.5' 1500/dx)\n"; 
  gp << "set xlabel 'x [km]'\n";
  gp << "set ytics out scale .5 rotate by 60 ('0.0' 0, '0.3' 300/dy, '0.6' 600/dy, '0.9' 900/dy, '1.2' 1200/dy, '1.5' 1500/dy)\n"; 
  gp << "set ylabel 'y [km]'\n";
  gp << "set output '" << file << "'\n";
  gp << "set grid\n";
  gp << "set multiplot layout " << ny << "," << nx << "\n";
}

template <class data_t>
void plot(Gnuplot &gp, const data_t &data)
{
  blitz::Array<float, 2> tmp(data);

  gp << "set xrange [0:" << tmp.extent(0)-1 << "]\n";
  gp << "set yrange [0:" << tmp.extent(1)-1 << "]\n";
  gp << "splot '-' binary" << gp.binfmt(tmp) << " origin=(0,0,0) with image failsafe notitle\n";

  gp.sendBinary(tmp);
}
