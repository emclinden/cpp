/* -*- C++ -*-
 *
 * ---------------------------------------------------------------------
 * $Id: applycalib.cpp 1062 2015-09-04 22:29:30Z mclinden $
 * ---------------------------------------------------------------------
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA
 *
 * ---------------------------------------------------------------------
 *
 */

/*

 *
 * Usage: applycalib result_file_from_fitcalibstars.dat
 
 */

//#define LTL_USE_SIMD


#include <ltl/marray.h>
#include <ltl/ascio.h>
#include <ltl/fitsio.h>
#include <ltl/misc/exceptions.h>
#include <ltl/util/command_line_reader.h>
#include <ltl/util/config_file_reader.h>
#include <ltl/util/option_parser.h>
#include <ltl/util/gnuplot.h>
#include <ltl/statistics.h>

#include <libcure/util.h>
#include <libcure/norm_g.h>

#include <libmath/interp.h>

using namespace util;

#include <sstream>
#include <cstdlib>
#include <string>

using std::string;
using ltl::MArray;
using ltl::FitsOut;
using std::exception;

const string version = "$Id: fitcalibstars.cpp 1062 2015-08-28 22:29:30Z mclinden $";

struct globals
{
      bool help;
      int sigma; //sigma value need to determine delta chi2
      string fitsfile;  //fits file
      string wave; //
} g;

void print_usage (OptionParser& flags);
void add_options (OptionParser& flags);



int main (int argc, char** argv)
{
   cout<<"[applycalib] VERSION: "<<version<<endl;

   // set options with commandline
   CommandLineReader comline (argc, argv); // read from command line
   OptionParser flags (&comline);
   add_options (flags);
   try
   {
      flags.parseOptions ();
   } catch (UException& e)
   {
      cerr<<"[fitcalibstars] ERROR: Caught exception: "<<e.what ()<<endl;
      print_usage (flags);
      exit (1);
   }

   if (g.help)
   {
      print_usage (flags);
      cout<<endl;
      return 0;
   }

   // main processing
   try
   {
      flags.parseOptions ();
      const list<string> files = comline.fileArguments (); // get the filenames
      if (files.size ()!=1)
         throw ltl::IOException ("[applycalib] ERROR: Exactly one results file from fitcalibstars must be given!");

      cout<<"[applycalib] MESSAGE: reading results from file '"<<(files.front ())<<"'"<<endl;

      ltl::AscFile in (files.front ());
      vector<string>modname;
      in.readColumn (1, modname); // model name
      MArray<float,1> c2 = in.readFloatColumn (2); //chi2 value
      MArray<float,1> ns = in.readFloatColumn (3); //number of samples
      float z;
      if (g.sigma == 1)
	{
	  z=0.68;
	}
      else if (g.sigma == 2)
	{
	  z = 1.39;
	}
      else if (g.sigma == 3)
	{
	  z = 2.09;
	}
      else
	{
	  throw std::invalid_argument ("expect sigma = 1, 2, or 3");
	}
      float dc2 = z * sqrt(2*ns(1)) * (c2(1)/ns(1));
      ltl::IndexList<1> w = where(c2 <= min(c2)+dc2);
      int n = w.size();
      MArray<float,1> modf;
      float c = 2.99792458e+18; //speed of light in A
      float lambda = 4694.35; //sdss g central wavelength
      float norm, fnu, mag, s, scale;
      MArray<float,2> ss (1024,n); //for scaled model spec (TODO: remove hardcoded 1024) 
      MArray<float,1> sc (n); //for scale arrays
      float endnum;

      ltl::FitsIn infits (g.fitsfile);
      MArray<float,2> fitsarray;
      infits>>fitsarray;

      ltl::AscFile win (g.wave);
      MArray<float,1> wvl;
      wvl = win.readFloatColumn (1);

      for (int i=1; i<=w.size(); ++i)
	{
	  modf = fitsarray(ltl::Range::all(),i);
	  norm = norm_to_sdss_g(wvl,modf);
	  fnu = (norm * lambda * lambda) / (c);
	  mag = -2.5*log10(fnu) - 48.6;
	  s = pow(10,((20.0 - mag)/2.5)); //TODO: remove hardcoded 20.0
	  scale = 1.0 / s;
	  ss(ltl::Range::all(),i)=modf*scale;
	  sc(i)=scale;

	}

     ltl::FitsOut outt("scaledmod_dchi2.fits");
     outt<<ss;    //write scaled model to fits
     ltl::FitsOut outf("scale_dchi2.fits");
     outf<<sc;    //write scale for each model to fits

   } catch (exception& e)
   {
      cerr<<"[applycalib] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
   return 0;
}


void print_usage (OptionParser& flags)
{
   cout<<"Usage: applycalib [options] results_file_from_fitcalibstars.cpp "<<endl;
   cout<<"  Use best fit models found by fitcalibstars to create sensitivity function."<<endl;
   cout<<endl;
   flags.printUsage (cout);
   exit (1);
}

void add_options (OptionParser& flags)
{
   try
   {
     flags.addOption (new IntOption ("sigma", "1", "include models within sigma of smallest chi2 (1, 2, or 3)", 's', &g.sigma));

      flags.addOption (new StringOption ("fitsfile", "models.fits", "fits file", 'f', &g.fitsfile));

      flags.addOption (new StringOption ("wave", "wvl.dat", "wavelength file", 'w', &g.wave));

   } catch (UException& e)
   {
      cerr<<"[applycalib] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
}
