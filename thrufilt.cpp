/* -*- C++ -*-
 *
 * ---------------------------------------------------------------------
 * mclinden 
 * ---------------------------------------------------------------------
 *
 * ---------------------------------------------------------------------
 *
 */

/*
 * get u-g colors for every model in fits file
 * remember: vectors[] are zero indexed, MArrays() are indexed at 1
 */
#include <ltl/marray.h>
#include <libmath/interp.h>
#include <ltl/ascio.h>
#include <ltl/fitsio.h>
#include <ltl/misc/exceptions.h>
#include <ltl/util/command_line_reader.h>
#include <ltl/util/config_file_reader.h>
#include <ltl/util/option_parser.h>
#include <string>

using ltl::MArray;
using ltl::FitsIn;
using namespace util;


struct globals
{
  string modlist;     //list of models that correspond to given fits file
  string wave;        //file containing wavelength array for templates
} g;

float thru_sdss_g (const MArray<float,1>& l, const MArray<float,1>& f);
float thru_sdss_u (const MArray<float,1>& l, const MArray<float,1>& f);
void print_usage (OptionParser& flags);
void add_options (OptionParser& flags);


float thru_sdss_g (const MArray<float,1>& l, const MArray<float,1>& f)
{
  //g filter wavelengths and transmission
   MArray<float,1> lg (18), g (18);
   lg = 3705.0f, 3830.0f, 3955.0f, 4080.0f, 4205.0f, 4330.0f, 4455.0f, 4580.0f, 4705.0f, 4830.0f, 4955.0f, 5080.0f, 5205.0f, 5330.0f, 5455.0f, 5580.0f, 5705.0f, 5830.0f;
   g = 0.0f, 0.0f, 0.1324f, 0.2609f, 0.3186f, 0.3496f, 0.3709f, 0.3863f, 0.3986f, 0.4043f, 0.4129f, 0.4201f, 0.4115f, 0.2112f, 0.0190f, 0.0024f, 0.0f, 0.0f;

   const int ng = g.nelements();
   const int xmin = find_closest_index (l, lg (1));
   const int xmax = find_closest_index (l, lg (ng));
   float sumf = 0, sumg = 0;
   for (int i = xmin; i<=xmax; ++i)
   {
      const float dw = (i>xmin+1 ? (l (i)-l (i-1)) : l (xmin+1)-l (xmin));
      float gint = lin_interp (g, lg, l (i)); //interpolate transmission to spectrum's wavelength array
      sumg += gint*dw;  
      sumf += gint*f (i)*dw;
   }
   return sumf/sumg; //approximate integration of flux in bandpass

}

float thru_sdss_u (const MArray<float,1>& l, const MArray<float,1>& f)
{
  //u filter wavelengths and transmission
  MArray<float,1> lu (15), u (15);
  lu = 2980.0f, 3005.0f, 3030.0f, 3155.0f, 3280.0f, 3405.0f, 3530.0f, 3655.0f, 3780.0f, 3905.0f, 4030.0f, 4055.0f, 4080.0f, 4105.0f, 4130.0f;
  u = 0.00f, 0.00140f, 0.00710f,  0.06290f,   0.13520f,   0.17370f,   0.1870f,   0.18410f,   0.14190f,  0.04050f, 0.00150f, 0.0008f, 0.0006f, 0.0003f, 0.000f;

   const int nu = u.nelements();
   const int xmin = find_closest_index (l, lu (1));
   const int xmax = find_closest_index (l, lu (nu));
   float sumf = 0, sumu = 0;
   for (int i = xmin; i<=xmax; ++i)
   {
      const float dw = (i>xmin+1 ? (l (i)-l (i-1)) : l (xmin+1)-l (xmin));
      float uint = lin_interp (u, lu, l (i)); //interpolate transmission to spectrum's wavelength array
      sumu += uint*dw;
      sumf += uint*f (i)*dw;
   }
   return sumf/sumu; //approximate integration of flux in bandpass

}

int main (int argc, char** argv)
{
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

   try
   {
   flags.parseOptions ();
   const list<string> file = comline.fileArguments (); // get the fits filename
   if (file.size ()!=1)
         throw ltl::IOException ("[thrufilt] ERROR: Only one fits file containing models can be given!");

   cout<<"[thrufilt] MESSAGE: reading list of models from "<<g.modlist<<endl;
   ltl::AscFile in(g.modlist);
   vector<string>modnames;
   in.readColumn (1,modnames); //read in (vector) of file names
   
   cout << "[thrufilt] MESSAGE: reading fits file '"<<(file.front())<<"'"<<endl;
   ltl::FitsIn infits(file.front());
   MArray<float,2> fitsarray;
   infits >> fitsarray; //read fits file and assign it to fitsarray

   ltl::AscFile win(g.wave);
   MArray<float,1> wvl;
   wvl = win.readFloatColumn (1); //native wavlength arrry of model
   float c = 2.99792458e+18; //speed of light in A
   float lambda_g = 4694.35; //sdss g central wavelength   
   float lambda_u = 3551.0; //sdss u central wavelength
   float nu, fnu_u, mag_u, ng, fnu_g, mag_g;
   string fout = "T03500_T10500_ug.dat"; //output file will contain model name and u-g color
   std::ofstream outf (fout.c_str ());  //TODO: change so fout is not hardcoded but depends on modlist
   outf << "#modname norm_u fnu_u mag_u norm_g fnu_g mag_g mag_u-mag_g"<<endl;

   //for each model calculate u - g and save it to file
   for (int i = 1; i<=modnames.size(); ++i) 
     {
       nu = thru_sdss_u(wvl,fitsarray(ltl::Range::all (),i)); //get f_lambda thru u filter
       fnu_u = (nu *lambda_u * lambda_u) / c; //convert f_lambda to f_nu
       mag_u = -2.5*log10(fnu_u); //convert f_nu to u magnitude
       ng = thru_sdss_g(wvl,fitsarray(ltl::Range::all (),i)); //get f_lambda thru g filter
       fnu_g = (ng *lambda_g * lambda_g) / c; //convert f_lambda to f_nu
       mag_g = -2.5*log10(fnu_g); //convert f_nu to g magnitude
       outf << modnames[i-1] << " " << mag_u - mag_g <<endl; //write to file
     }
   

   } catch (exception& e)
   {
      cerr<<"[fitcalibstars] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
   return 0;
}


void print_usage (OptionParser& flags)
{
  cout<<"Usage: thrufilt -m model_list.txt fitsfile_models.fits"<<endl;
   cout<<endl;
   flags.printUsage (cout);
   exit (1);
}

void add_options (OptionParser& flags)
{
   try
   {
      flags.addOption (new StringOption ("modlist", "T03500_T10500.list", "Input file containing model names corresponding to supplied fitsfile.", 'm', &g.modlist));

      flags.addOption (new StringOption ("wave", "LAMBDA_D1A_long.dat", "Input file containing wavelength array for templates.", 'w', &g.wave));

   } catch (UException& e)
   {
      cerr<<"[thrufilt] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
}


