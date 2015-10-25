/* -*- C++ -*-
 *
 * ---------------------------------------------------------------------
 * $Id: fitcalibstars.cpp 1062 2015-08-28 22:29:30Z mclinden $
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
 * Usage: fitcalibstars [options] -d ditherfile.txt -i IFUfile -t template_list detection_cont.txt
 
 */

//#define LTL_USE_SIMD


#include <ltl/marray.h>
#include <ltl/ascio.h>
#include <ltl/fitsio.h>
#include <ltl/misc/exceptions.h>
#include <libcure/util.h>

#include <libcure/ditherenv.h>
#include <libcure/extractor_prof.h>

#include <ltl/util/command_line_reader.h>
#include <ltl/util/config_file_reader.h>
#include <ltl/util/option_parser.h>

#include <ltl/statistics.h>
#include <ltl/convolve.h>
#include <ltl/util/gnuplot.h>
#include <ltl/util/timer.h>

#include <libmath/bspline.h>
#include <libmath/interp.h>

using namespace util;

#include <sstream>
#include <cstdlib>
#include <string>

using std::string;
using ltl::FitsOut;
using ltl::FitsIn;
using ltl::MArray;

const string version = "$Id: fitcalibstars.cpp 1062 2015-08-28 22:29:30Z mclinden $";

struct globals
{
  bool help;
  string outfile;     // output image name
  string prefix;      // output image prefix
  string dither_file; // dither description
  string ifu_file;    // IFU description
  int nsamples;       // number of wavelength samples
  float lambda[2];    // wavelength range
  bool fx;            // do not resample wavelengths
  bool plot;          //plot results?
  float kron;         // aperture scaling factor factor
  float minap;        // Minimal aperture size
  float aperture;     // extraction aperture
  string templates;   // file containing list of Castelli+Kurucz templates 
  string colorfile;   //file containing u-g color for templates
  float color;        //observed u - g color
  float ecolor;       //error on observed u - g color
  float snr;          // snr to consider fiber
  string wave;        //file containing wavelength array for templates
  string fitsfile;    // fits file containing model templates
  bool fitswrite;     //whether to write out fits files containing resampled data and splines
  bool txtwrite;      //whether to write out txt files containing resampled data and splines
  bool which;         //use starextract output if true
  string starex;      //output from starextract needed if which true
  string writeprefix; //prefix for fits files if fitswrite or txtwrite on
  int nresults;       //number of results files to write out
} g;

void print_usage (OptionParser& flags);
void add_options (OptionParser& flags);
float int_trapezium (const MArray<float,1>& Spec, float from, float to, float n);
MArray<float,1> medfilt(const MArray<float,1>& inspec, int window);

struct MySort
{
      MArray<float,1> myarray;

      MySort (MArray<float,1> & val_vec)
      {
         myarray = val_vec;
      }

      bool operator() (int i1, int i2)
      {
         return myarray (i1)<myarray (i2);
      }
};

struct MyResults
{
      MArray<float,1> chi2;
      vector<string> tname;
};

struct GaussianKernel1D
{
      GaussianKernel1D (const float s, const float k) :
            sigma_ (s), extent_ (k)
      {
         g_.makeReference (MArray<float,1> (2*extent_+1));
         g_.setBase (-extent_);
         g_ = exp (-0.5*(pow2 (indexPosDbl (g_, 1)))/(sigma_*sigma_));
         g_ /= sum (g_);
      }

      template<typename Iter>
      inline typename Iter::value_type eval (Iter& A) const
      {
         float r = 0.0;
         for (int i = -extent_; i<=extent_; ++i)
            r += g_ (i)*OFFSET1(A, i);

         return typename Iter::value_type (r);
      }

      float sigma_;
      int extent_;
      MArray<float,1> g_;
};

namespace ltl
{
   template<typename Iter>
   struct kernel_return_type<Iter,GaussianKernel1D>
   {
      typedef typename Iter::value_type value_type;
   };
}

float int_trapezium (const MArray<float,1>& Spec, float from, float to, float n)
{
   float h = (to-from)/n;
   int i;
   float sum = (h/2.0)*(Spec (1)+Spec (n));
   for (i = 2; i<n; i++)
      sum += h*Spec (i);
   return sum;
}

MArray<float,1> medfilt(const MArray<float,1>& inspec, int window)
{
  int n = inspec.size();
  int addelems = (window+1)/2;
  MArray<float,1> newarray (n+(2*addelems));

  //add buffers to ends of array, size of buffer on each end is windowsize/2
  for (int i = 1; i <= addelems; ++i)
    {
      newarray(addelems + 1 - i)=inspec(i);
      newarray(n + addelems +i ) = inspec(n + 1 - i);
    }

  newarray(ltl::Range(addelems+1,n+addelems))=inspec; //fill in the middle of newarray with inspec

  MArray<float,1> outarray (n+2); //empty array to hold output
  float temp = 0.0;
  for (int i = 1; i <= newarray.size() - window; ++i)
    {
      temp = median_exact(newarray(ltl::Range(i,i+window-1))); //what's ltl marray for median?
      outarray(i) = temp;
    }

   MArray<float,1> finalarr (n);
   //slice off first and last entries in outarray
   finalarr(ltl::Range::all())=outarray(ltl::Range(2,outarray.size()-1));
   return finalarr;
}

std::string tostring (int toconv)
{
   stringstream ss;
   ss<<toconv;
   string toret (ss.str ());
   return toret;
}

MArray<float,2> extract_fibers_in_aperture (const std::vector<ResElemList>& data, const DitherEnvironment& DE, int id)
{
   MArray<float,1> wave (g.nsamples);
   int counter = 0;
   MArray<float,1> totspec (g.nsamples);
   totspec (ltl::Range::all ()) = 0.0;
   MArray<float,1> totespec (g.nsamples);
   totespec (ltl::Range::all ()) = 0.0;
   // loop over all frames/dithers/shots
   for (unsigned int s = 0; s<data.size (); ++s)
   {
      const ResElemList& shot = data[s];
      // loop over all resolution elements
      for (ResElemList::const_iterator p = shot.begin (); p!=shot.end (); ++p)
      {
         const ResElem& res = *p;
         if (!res.ignore)
         {

            const MArray<float,2> & Data = DE.get_vframe (s).get_data (res.tu); //get data
            const MArray<float,2> & Err = DE.get_vframe (s).get_edata (res.tu); //get errors
            const Distortion& D = DE.get_vframe (s).get_distortion (res.tu);
            const float lam1 = D.map_xy_wavelength (__minx__+5, (__maxy__-__miny__)/2);
            const float lam2 = D.map_xy_wavelength (__maxx__-5, (__maxy__-__miny__)/2);
            MArray<float,2> Err2;
            Err2 = Err;
            Err2 = pow2 (Err2);

            Extractor E (D);

            MArray<float,2> S2 (g.nsamples, 1); //to hold data
            S2 = 0.0f;
            MArray<float,2> E2 (g.nsamples, 1); //to hold errors
            E2 = 0.0f;
            //for a given shot and resolution element, extract data and errors in the relevant fiber
            S2 (ltl::Range::all (), 1) = E.extract_aperture (Data, res.tf, g.aperture, lam1, lam2, g.nsamples);
            E2 (ltl::Range::all (), 1) = E.extract_aperture (Err2, res.tf, g.aperture, lam1, lam2, g.nsamples);
            E2 = sqrt (E2);

            //create some strings to declare which shot, ifu, and fiber you're looking at
            std::string outfib = tostring (res.tf);
            int shotout = s+1;
            string outifu = (res.tu+1==1 ? "L" : "R");
            std::string outshot = tostring (shotout);
            std::string outobj = tostring (id);
            string outname;
            string eoutname;
            outname = "obj"+outobj+"_D"+outshot+"_"+outifu+"_fiber_"+outfib+".fits"; //output fits file name
            eoutname = "e."+outname;

            FitsOut O (outname);
            O.addValueCard ("CTYPE1  ", "LINEAR", ""); //add wavelength info
            O.addValueCard ("CRPIX1  ", 1.0f, "");
            O.addValueCard ("CUNIT1  ", "Angstrom", "");
            O.addValueCard ("CRVAL1  ", lam1, "");
            O.addValueCard ("CDELT1  ", (lam2-lam1)/(g.nsamples-1), "");
            O.addValueCard ("CD1_1   ", (lam2-lam1)/(g.nsamples-1), "");
            O<<S2;  //write data to fits file

            FitsOut EO (eoutname);
            O.addValueCard ("CTYPE1  ", "LINEAR", ""); //add wavelength info
            O.addValueCard ("CRPIX1  ", 1.0f, "");
            O.addValueCard ("CUNIT1  ", "Angstrom", "");
            O.addValueCard ("CRVAL1  ", lam1, "");
            O.addValueCard ("CDELT1  ", (lam2-lam1)/(g.nsamples-1), "");
            O.addValueCard ("CD1_1   ", (lam2-lam1)/(g.nsamples-1), "");
            EO<<E2;  //write error to fits file
 
            MArray<float,1> a (g.nsamples);
            MArray<float,1> ea (g.nsamples);
	    a = S2 (Range::all (), 1);
	    ea = E2 (Range::all (), 1);

            const float cd = (lam2-lam1)/(g.nsamples-1);
            MArray<float,1> w (g.nsamples); //
            w = (indexPosDbl (w, 1)-1.0)*cd+lam1;

            float trap = int_trapezium (a, w (1), w (g.nsamples), w.size ());
            float etrap = int_trapezium (ea, w (1), w (g.nsamples), w.size ());
            float snr = trap/etrap;
            if (snr>=g.snr)
            {
               counter++;
               totspec += a;
               totespec += ea*ea;
            }
            wave (ltl::Range::all ()) = w;
         }
      }
   }

   MArray<float,2> ret (g.nsamples, 3);
   ret (ltl::Range::all (), 1) = totspec/counter;
   ret (ltl::Range::all (), 2) = sqrt (totespec)/counter;
   ret (ltl::Range::all (), 3) = wave;
   return ret;
}

MyResults find_best_mod (MArray<float,2>& tofit, int k, int id)
{
   MArray<float,1> a (g.nsamples);   // obs spec before filtering
   MArray<float,1> aa (g.nsamples);   // obs spec after filtering
   MArray<float,1> ea (g.nsamples);  // errors
   MArray<float,1> w (g.nsamples);   // wavelength grid
   MArray<float,1> sa (g.nsamples);  // obs continuum spline evaluated at w
   MArray<float,1> ma (g.nsamples);  //median filtered obs spec

   a = tofit (Range::all (), 1);
   ea = tofit (Range::all (), 2);
   w = tofit (Range::all (), 3);
   
   ma = medfilt(a,9);  //median filter obs spec to find cosmics (use window=9 for pickles, window=3 for kurucz)
   aa = merge(a-ma > 1000, ma, a);  //where difference is large, replace with medianed value (120)
   ltl::IndexList<1> c2bad = where(a-ma > 1000); //keep track of where spec altered to remove during chi2 calculation (120)

   // spline fit to continuum of observed spectrum
   BSpline::BSpline<float> S (w.beginRA (), w.nelements (), aa.beginRA (), 1900, BSpline::BSplineBase<float>::BC_ZERO_SECOND);

   //divide averaged fiber spectrum by bspline curve
   MArray<float,1> anorm (g.nsamples);
   MArray<float,1> danorm (g.nsamples);  //TODO: avoid Inf?

   // divide observed spectrum by spline and set errors to large value 
   // instead of masked values
   sa = apply(S,w);
   anorm = aa/sa;   
   danorm = merge (ea>0.0, ea/sa, 1e6);

   //read wavelength array for templates
   ltl::AscFile win (g.wave);
   MArray<float,1> wx;
   wx = win.readFloatColumn (1); //native wavlength array of model
   //read in list of model files
   vector<string>filename;
   ltl::AscFile in (g.templates);
   in.readColumn (1, filename);
   int sz = filename.size ();

   ltl::AscFile uin (g.colorfile); //read color info for templates
   MArray<float,1> ug;
   ug = uin.readFloatColumn (2);
   
   //set up empty arrays for model spectra
   MArray<float,1> y;
   MArray<float,1> ylin;
   GaussianKernel1D G (7.44,21); //create gaussian kernel for convolution
   MArray<float,1> reschi2 (sz);
   MArray<float,2> mods (g.nsamples, sz);
   MArray<float,2> modsp (g.nsamples, sz);
   MArray<float,1> denom (sz); 
   MArray<float,1> resalpha (sz);

   util::Timer t1;
   t1.start();
   ltl::FitsIn infits (g.fitsfile);
   MArray<float,2> fitsarray;
   infits>>fitsarray;
   t1.stop ();
   cout<<"read fits "<<t1.elapsedSeconds () << endl;

   t1.start ();

   //read spectrum of each model
   for (int i = 1; i<=filename.size (); ++i)
   {
      y = fitsarray (ltl::Range::all (), i);
      MArray<float,1> C (wx.size ());
      C = 0.0f;
      C = convolve (y, G); //convolve y array with gaussian
      ylin = lin_interp_ff (C, wx, w); //linearly interpolate to same wavelength grid

      //calculate bspline for continuum shape of model
      BSpline::BSpline<float> M (w.beginRA (), w.nelements (), ylin.beginRA (), 1900, BSpline::BSplineBase<float>::BC_ZERO_SECOND);

      //divide model spec by continuum spline
      MArray<float,1> mnorm (g.nsamples);
      MArray<float,1> meval (g.nsamples);
      meval = apply(M, w);
      mnorm = ylin/meval;

      float alpha = int_trapezium(anorm,w(1),w(g.nsamples),w.size()) / int_trapezium(mnorm,w(1),w(g.nsamples),w.size());
      
      //calculate x^2 of this model
      MArray<float,1> c2 (g.nsamples); 
      c2 = pow2(anorm - (alpha*mnorm)) / ((danorm*danorm) + (0.0016)); 
      if (c2bad.size() != 0) //if there were cosmics removed, don't use these pixels to calculate c2
	{
	  c2(c2bad)=0;
      	}     
  
      float chi2m = sum(c2,0.0);
      reschi2(i)=chi2m/(g.nsamples - c2bad.size() - 1);
      denom(i)=g.nsamples - c2bad.size () -1; 
      mods (ltl::Range::all (), i) = ylin; //conv, resamp model
      modsp (ltl::Range::all (), i) = meval; //model spline
      resalpha(i)=alpha;
      
   }
   t1.stop ();
   cout<<"fit models "<<t1.elapsedSeconds () << endl; //how long it took to fit models

   MArray<int,1> ww (sz);
   ww = indexPosInt (w, 1); //set up array of integer indices
   MArray<float,1> chi2_ret (sz); //for sorted chi2
   vector<string> tname_ret (sz); //for sorted names
   MArray<float,2> mods_ret (g.nsamples, sz); //for sorted model spec
   MArray<float,2> modsp_ret (g.nsamples, sz); //for sorted model splines

   //sort chi2 array from min to max - use this to sort indices of filename vector, etc
   std::sort (ww.beginRA (), ww.endRA (), MySort (reschi2)); //sort chi2
   string rout = g.dither_file+"results_obj"+tostring (id)+"_"+g.templates+".txt";
   std::ofstream outr (rout.c_str ());

   for (int i = 1; i<=sz; ++i)
   {
      chi2_ret (i) = reschi2 (ww (i));
      tname_ret[i-1] = filename[ww (i)-1];
      string yesno;
      (ug (ww(i)) >= g.color - g.ecolor and ug(ww(i)) <= g.color + g.ecolor) ? yesno = "yes" : yesno = "no";
      //save the sorted chi2 and filename arrays to results file
      outr<<filename[ww (i)-1]<<" "<<reschi2 (ww (i))<<" "<< denom (ww (i))<<" "<<resalpha(ww(i))<<" "<<yesno<<endl; 
      //sort models, splines, input spec etc based on chi2
      mods_ret (ltl::Range::all (), i) = mods (ltl::Range::all (), ww (i));
      modsp_ret (ltl::Range::all (), i) = modsp (ltl::Range::all (), ww (i));
   }

   //if writeout on, return (big) fits files for models, model splines, and observed spectra
   //this is slow - only turn on if really needed
   if (g.fitswrite)
   {
     ltl::FitsOut outm(g.writeprefix+"_models.fits");
     ltl::FitsOut outs(g.writeprefix+"_mspline.fits");
     ltl::FitsOut outo(g.writeprefix+"_observed.fits");
     MArray<float,2> mnout(1024,67838);
     MArray<float,2> msout(1024,67838);
     MArray<float,2> oout(1024,5);
     for (int i = 1; i<=filename.size(); ++i)
      {
     	mnout(ltl::Range::all(),i) = mods_ret(ltl::Range::all(),i);
     	msout(ltl::Range::all(),i) = modsp_ret(ltl::Range::all(),i);
      }
     outm << mnout;
     outs << msout;
     oout(ltl::Range::all(),1)=w;
     oout(ltl::Range::all(),2)=a;   
     oout(ltl::Range::all(),3)=aa;
     oout(ltl::Range::all(),4)=ea;
     oout(ltl::Range::all(),5)=sa;
     outo << oout;
   }
   //turn on if you want ascii files returned with spectra, splines, errors
   if (g.txtwrite)
   {
     for (int i = 1; i<=g.nresults; ++i)
      {
         string pout = "obj"+tostring (id)+"_"+tname_ret[i-1]+".out";
         std::ofstream outp (pout.c_str ());
	 outp<<"#wavelength "<<"model "<<"model_spline "<<"filtered_obs_spectrum "<<"unfiltered_obs_spectrum_"<<"obs_error "<<"obs_spline "<<endl;
         for (int j = 1; j<=g.nsamples; ++j)
         {
	   outp<<w (j)<<" "<<mods_ret (j, i)<<" "<<modsp_ret (j, i)<<" "<<aa(j)<<" "<<a(j)<<" "<<ea(j)<<" " <<sa (j)<<endl;
         }
      }
   }

   //return sorted chi2 and sorted names of model
   struct MyResults returnthis;
   returnthis.chi2 = chi2_ret (ltl::Range::all ());
   returnthis.tname = tname_ret;
   return returnthis;
}

int main (int argc, char** argv)
{
   cout<<"[fitcalibstars] VERSION: "<<version<<endl;

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
         throw ltl::IOException ("[fitcalibstars] ERROR: Exactly one detection file must be given!");

      // read catalog data id, coordinates, and apertures
      cout<<"[fitcalibstars] MESSAGE: reading catalog from file '"<<(files.front ())<<"'"<<endl;
      ltl::AscFile in (files.front ());
      MArray<int,1> id = in.readIntColumn (1);
      MArray<float,1> x = in.readFloatColumn (2);
      MArray<float,1> y = in.readFloatColumn (3);
      MArray<float,1> a = in.readFloatColumn (6);
      MArray<float,1> b = in.readFloatColumn (7);
      MArray<float,1> pa = in.readFloatColumn (8);
      MArray<float,1> ka = in.readFloatColumn (10);
      MArray<float,1> kb = in.readFloatColumn (11);

      DetectOptions DO;
      DitherEnvironment DE;
      ProjectorSet PS;

      const int nobj = id.length (1);
      vector<string>sanames;

      if (g.which) 
	{
	  //if using star extracted files read  list fits files here
	  cout<<"[fitcalibstars] MESSAGE: reading list of starextract fits files from '"<<(g.starex)<<"'"<<endl;
	  ltl::AscFile in (g.starex);
	  in.readColumn (1, sanames);
	  if (nobj != sanames.size())
	    throw ltl::IOException ("[fitcalibstars] ERROR: Number of objects in detect file and number of files in starextract list must match!");
	}
      else
	{
	  //if not using starextract results, setup detectoptions, ditherenv, projectorset once outside of for-loop
	  DO.psf_size[0] = DO.psf_size[1] = 6;
	  DO.sub_cont = false;
	  DO.ifu_file = g.ifu_file;

	  DE.set_detect_parameters (DO);
	  DE.set_cmd_line (argc, argv);
	  DE.read_ditherfile (g.dither_file);
	  DE.setup_projectors (true);

	  PS = DE.get_projector_set ();
	}
      
      MArray<float,2> ans (g.nsamples,3);
      for (int i = 1; i<=nobj; ++i)
      {
	if (g.which)
	  {
	    //get stellar spectrum, error array and wvl from starextract output files
	    cout<<"[fitcalibstars] MESSAGE: reading output spectrum from starextract from file "<<sanames[i-1]<<endl;

	    ltl::FitsIn starfits (sanames[i-1]); 
	    MArray<float,2> starfitsarray;
	    starfits>>starfitsarray;   //stellar spectrum
	    MArray<float,1> sfa = starfitsarray(ltl::Range::all(),1);
	    sfa = merge(sfa != sfa, 0., sfa);
	    ltl::IndexList<1> fil = where(sfa != sfa);
	    sfa(fil) = 0; //where starextract returned nan, substitute 0

	     ltl::FitsIn estarfits ("e."+sanames[i-1]);
             MArray<float,2> estarfitsarray;
             estarfits>>estarfitsarray;  //error spectrum
	     MArray<float,1> esfa = estarfitsarray(ltl::Range::all(),1);
	     cout<<"[fitcalibstars] MESSAGE: reading error spectrum from "<<"e."+sanames[i-1]<<endl;

	     //create wvl array from header info
             float crv = starfits.getFloat ("CRVAL1  ");
             float cd1 = starfits.getFloat ("CD1_1   ");
             MArray<float,1> wvl(g.nsamples);
             wvl = (indexPosDbl (wvl, 1)-1.0)*cd1+crv; //wavelength array

	     //save everything in one structure called ans for find_best_mod
	     ans (ltl::Range::all (), 1) = sfa;
	     ans (ltl::Range::all (), 2) = esfa;
	     ans (ltl::Range::all (), 3) = wvl;	   
	  }
	else
	  {
	    //if not using starextract results, go through finding fibers
	    cout<<"[fitcalibstars] MESSAGE: Extracting source "<<i<<" at wavelength of 4500A"<<endl;
	    cout<<"[fitcalibstars] MESSAGE:    aperture "<<2*a (i)*g.kron<<endl;
	    cout<<"[fitcalibstars] MESSAGE: aperture scale factor "<<g.kron<<endl;
	    float w = 4500;
	    std::vector<ResElemList> modelPS;
	    ProjectorSet::SkyPosition p (x (i), y (i), w);
	    const float diam = (2*a (i)*g.kron>g.minap ? a (i) : g.minap);
	    PS.findResElemsForPointSourceInAperture (modelPS, p, diam*g.kron);
	    ans = extract_fibers_in_aperture (modelPS, DE, id (i));
	  }

         struct MyResults fanswer; //structure that houses the results (2 marrays containing chi2 values and filenames)
	 
	 //check if some fibers have been found with s/n greater than threshold
	 //if so, go on to find best model
	 float keepgoing = sum (ans (ltl::Range::all (), 1));
         if (keepgoing>0.0) 
         {
            fanswer = find_best_mod (ans, i, id (i));
         }
         else
            cout<<"No fibers found above your s/n threshold, lower your threshold (or maybe don't use this star for calib)."<<endl;

         if (g.plot==1) //this only works if templates names are ck-like with 'T' in the name
         {
            //get temps of ck models
            MArray<float,1> temps (fanswer.chi2.size ());
            std::size_t pos1 = g.templates.find ("."); //split template file list name on '.'
            std::string gt = g.templates.substr (0, pos1); //keep first part of split for output file name
            for (int j = 1; j<=fanswer.chi2.size (); ++j)
            {
               std::string s = fanswer.tname[j+1];
               std::size_t pos2 = s.find ("T"); //find 'T' in template name
               std::string ss = s.substr (pos2+1, 5);
               const char * c = ss.c_str ();
               temps (j) = std::atof (c); //read next 5 chars after 'T' as the model temperature

            }
            ltl::Gnuplot gp;
            gp<<"set terminal pdf"<<endl;
            gp<<"set output \""<<"obj"+tostring (id (i))+"_chi2r_temp_"<<gt<<".pdf\""<<endl;
            gp<<"set xlabel \"Temperature\""<<endl;
            gp<<"set ylabel \"reduced chi2\""<<endl;
            if (max (fanswer.chi2)>=1000.)
               gp<<"  set yrange [0:700]"<<endl;
            else
               gp<<"set yrange ["<<0<<":"<<max (fanswer.chi2)+50<<"]"<<endl;
            gp<<"set nokey"<<endl;
            gp<<"plot '-' with points\n";
            gp.send (temps, fanswer.chi2);

         }

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
   cout<<"Usage (using starextract output): fitcalibstars [options] -S list_of_starextractfits.txt -t list_of_models.txt -f fitsfile_containing_models.fits detection_cont.txt "<<endl;
   cout<<"or"<<endl;
   cout<<"Usage (fit fibers in situ): fitcalibstars [options] -X -d ditherfile.txt -i IFUfile -t list_of_models.txt -f fitsfile_containing_models.fits detection_cont.txt "<<endl;
   cout<<"   Find best fit models for calibration stars on IFU."<<endl;
   cout<<endl;
   flags.printUsage (cout);
   exit (1);
}

void add_options (OptionParser& flags)
{
   try
   {
      flags.addOption (new BoolOption ("help", "FALSE", "Print usage information.", 'h', &g.help));

      flags.addOption (new StringOption ("outfile", "", "Output image base filename [=input image]", 'o', &g.outfile));

      flags.addOption (new StringOption ("prefix", "FF", "Output image filename prefix.", 'p', &g.prefix));

      flags.addOption (
            new StringOption ("ditherfile", "dither.txt", "The ditherfile describing the dithered observation if !g.which.", 'd', &g.dither_file));

      flags.addOption (new StringOption ("IFU-file", "IFUcen.txt", "Input file with fiber positions if !g.which.", 'i', &g.ifu_file));

      flags.addOption (new FloatOption ("ap-factor", "1.0", "Rescale default aperture size (FWHM from ditherfile).", 'k', &g.kron));

      flags.addOption (new FloatOption ("minimum-aperture", "3.0", "Minimal aperture size in arcsec.", 'm', &g.minap));

      flags.addOption (new BoolOption ("x", "FALSE", "Do not resample in wavelength.", 'x', &g.fx));

      flags.addOption (new BoolOption ("plot", "FALSE", "plot results?", 'P', &g.plot));

      flags.addOption (new IntOption ("nsamples", "1024", "Number of wavelength samples.", 'n', &g.nsamples));

      flags.addOption (new FloatArrayOption ("lambda", "0,0", "Wavelength range ([0,0] for full range).", 'l', 2, g.lambda));

      flags.addOption (new FloatOption ("aperture", "5.0", "Extraction aperture.", 'a', &g.aperture));

      flags.addOption (new FloatOption ("snr", "10.0", "S/N threshold to include fiber.", 's', &g.snr));

      flags.addOption (new BoolOption ("which", "TRUE", "use starextract results instead of finding and extracting high s/n fibers", 'X', &g.which));

      flags.addOption (new StringOption ("starextract-file", "starextract.list", "if g.which true need to supply this ascii file listing fits files from starextract", 'S', &g.starex));

      flags.addOption (new StringOption ("wave", "LAMBDA_D1A_short.dat", "Input file containing wavelength array for templates.", 'w', &g.wave));
      
      flags.addOption (new StringOption ("colorfile", "T03500_T10500_ug.dat", "Input file containing u-g colors for templates.", 'u', &g.colorfile));

      flags.addOption (new FloatOption ("color", "0.0", "u - g color of observed calibration star", 'c', &g.color));

      flags.addOption (new FloatOption ("ecolor", "0.0", "error onu - g color of observed calibration star", 'e', &g.ecolor));

      flags.addOption (new StringOption ("fitsfile", "T03000_T10000.fits", "Input fits file containing model templates.", 'f', &g.fitsfile));

      flags.addOption ( new BoolOption ("fitswrite", "FALSE", "write out (BIG!) fits file containing resampled models, data & splines.", 'F', &g.fitswrite));
      
      flags.addOption ( new BoolOption ("txtwrite", "FALSE", "write out ascii files containing resampled models, data & splines.", 'W', &g.txtwrite));

      flags.addOption (new StringOption ("writeprefix", "out", "prefix for fits files if either fitswrite or txtwrite on", 'R', &g.writeprefix));

      flags.addOption (new IntOption ("nresults", "100", "If txtwrite is on, how many results do you want files for?", 'r', &g.nresults));


      flags.addOption (
            new StringOption (
                  "templates", "templates.txt",
                  "list of Castelli+Kurucz templates names", 't',
                  &g.templates));

   } catch (UException& e)
   {
      cerr<<"[fitcalibstars] ERROR: caught exception: "<<e.what ()<<endl;
      exit (1);
   }
}
