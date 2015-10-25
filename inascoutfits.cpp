/* -*- C++ -*-
 *
 * ---------------------------------------------------------------------
 * $Id:$ inascoutfits.cpp mclinden $
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
#define LTL_RANGE_CHECKING

#include <ltl/marray.h>
#include <ltl/ascio.h> 
#include <ltl/misc/exceptions.h>
#include <numeric>
#include <ltl/fitsio.h>
#include <ltl/marray_io.h>

#include <libcure/util.h>

#include <ltl/util/command_line_reader.h>
#include <ltl/util/config_file_reader.h>
#include <ltl/util/option_parser.h>
#include <dirent.h>
#include <ltl/statistics.h>
#include <ltl/marray/reductions.h>


#include <sstream>
#include <cstdlib>
#include <string>
#include <list>

using std::string;
using std::list;
using ltl::FitsOut;
using ltl::FitsIn;
using ltl::MArray;
using ltl::Range;
using ltl::sum;
using namespace util;
using namespace ltl;

struct globals
{
  string outfile;     // output fits name
} g;



int main( int argc, char** argv )
{
  CommandLineReader comline(argc, argv); // read from command line

   try
   {
     size_t start = time(NULL);

     ltl::FitsIn in3("T03.list.fits");
     MArray<float,2> three;
     in3>>three;
     ltl::FitsIn in4("T04.list.fits");
     MArray<float,2> four;
     in4>>four;
     ltl::FitsIn in5("T05.list.fits");
     MArray<float,2> five;
     in5>>five;
     ltl::FitsIn in6("T06.list.fits");
     MArray<float,2> six;
     in6>>six;
     ltl::FitsIn in7("T07.list.fits");
     MArray<float,2> seven;
     in7>>seven;
     ltl::FitsIn in8("T08.list.fits");
     MArray<float,2> eight;
     in8>>eight;
     ltl::FitsIn in9("T09.list.fits");
     MArray<float,2> nine;
     in9>>nine;
     ltl::FitsIn in10("T10.list.fits");
     MArray<float,2> ten;
     in10>>ten;

     int aa =in3.getNaxis(2);
     int bb=aa+in4.getNaxis(2);
     int cc=bb+in5.getNaxis(2);
     int dd=cc+in6.getNaxis(2);
     int ee=dd+in7.getNaxis(2);
     int ff=ee+in8.getNaxis(2);
     int gg=ff+in9.getNaxis(2);
     int hh=gg+in10.getNaxis(2);

     ltl::FitsOut outt("T03000_T10000_v3.fits");
     MArray<float,2> fout(2201,hh);

     //shorten wavelength range (indices 901-3101) to make output file much smaller
     fout(ltl::Range::all(),ltl::Range(1,aa)) = three(ltl::Range(901,3101),ltl::Range::all());
     cout<<"three"<<endl;
     fout(ltl::Range::all(),ltl::Range(aa+1,bb)) = four(ltl::Range(901,3101),ltl::Range::all());
     cout<<"four"<<endl;
     fout(ltl::Range::all(),ltl::Range(bb+1,cc)) = five(ltl::Range(901,3101),ltl::Range::all());
     cout<<"five"<<endl;
     fout(ltl::Range::all(),ltl::Range(cc+1,dd)) = six(ltl::Range(901,3101),ltl::Range::all());
     cout<<"six"<<endl;
     fout(ltl::Range::all(),ltl::Range(dd+1,ee)) = seven(ltl::Range(901,3101),ltl::Range::all());
     cout<<"seven"<<endl;
     fout(ltl::Range::all(),ltl::Range(ee+1,ff)) = eight(ltl::Range(901,3101),ltl::Range::all());
     cout<<"eight"<<endl;
     fout(ltl::Range::all(),ltl::Range(ff+1,gg)) = nine(ltl::Range(901,3101),ltl::Range::all());
     cout<<"nine"<<endl;
     fout(ltl::Range::all(),ltl::Range(gg+1,hh)) = ten(ltl::Range(901,3101),ltl::Range::all());
     cout<<"ten"<<endl;
     outt <<fout;
     
     size_t end = time(NULL);
     double elapsed = double(end -start);
     cout<<elapsed<<endl;

   } catch (exception& e)
   {
      cerr<<"[findfiber] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
   return 0;
}

