/***************************************************************************
 *
 *   Copyright (C) 2010 by Paul Demorest
 *   Licensed under the Academic Free License version 2.1
 *
 ***************************************************************************/

using namespace std;

#include "Pulsar/Application.h"
#include "Pulsar/StandardOptions.h"
#include "Pulsar/Archive.h"
#include "Pulsar/Integration.h"
#include "Pulsar/Profile.h"
#include "Pulsar/DynamicSpectrum.h"
#include "Pulsar/StandardFlux.h"
#include "Pulsar/ProfileAmps.h" //added this

using namespace Pulsar;
using namespace std;

//
//! Produce dynamic spectra from Archives
//
class psrflux : public Pulsar::Application
{
public:

  //! Default constructor
  psrflux ();

  //! One-time setup
  void setup ();

  //! Process the given archive
  void process (Pulsar::Archive*);

  //! Extension to append
  std::string ext;

  //! Flux calculation
  //DynamicSpectrum ds;

  vector <DynamicSpectrum> ds;

  //! For profile aligment if needed
  ProfileShiftFit psf;

  //! Standard profile file
  std::string stdfile;

  //! Standard archive
  Reference::To<Archive> stdarch;

  //! Setup to use the given file as a standard
  void set_standard(Archive *arch);

  //! Align w/ each input profile separately
  bool align;

  //! Disable fits for alignment
  bool noalign;

  //compute lcp and rcp dynamic spectra 0 and 1 
  bool separate_cp_spectra;

protected:

  //! Add command line options
  void add_options (CommandLine::Menu&);

  //! Standard preprocessing options
  Pulsar::StandardOptions standard_options;

};


psrflux::psrflux ()
  : Application ("psrflux", "Produce dynamic spectra from Archives")
{
  ext = ".dynspec";
  stdfile = "";
  align = false;
  noalign = false;

  ds.resize(1);

  add( &standard_options );
}

void psrflux::add_options (CommandLine::Menu& menu)
{
  CommandLine::Argument* arg;

  menu.add("\n" "Dynamic spectrum options. You've been changed!:");

  arg = menu.add (ext, 'e', "ext");
  arg->set_help ("Append extention to output (default .ds)");

  arg = menu.add (stdfile, 's', "std");
  arg->set_help ("Standard profile file");

  arg = menu.add (align, 'a', "align");
  arg->set_help ("Align standard with each profile separately");

  arg = menu.add (noalign, 'A', "noalign");
  arg->set_help ("No fit for profile alignment at all");

  arg = menu.add (separate_cp_spectra, 'c');
  arg-> set_help("Compute LCP RCP dynspec");
}

void psrflux::set_standard(Archive *arch)
{
  if(separate_cp_spectra){
  cout << "Setting standard" << endl;
  ds.resize(2);

  stdarch = arch -> clone();
  stdarch -> tscrunch();
  stdarch -> fscrunch();

  const unsigned nbin = stdarch -> get_nbin();
  const unsigned isub = stdarch -> get_nsubint();
  cout << "nbin: " << nbin << endl;
  cout << "nsubint: " << isub << endl;
  
  stdarch -> convert_state(Signal::Stokes);

  Pulsar::Integration* integration = stdarch -> get_Integration(0);  
  
  Pulsar::Profile* I = integration -> get_Profile(0,0);
  Pulsar::Profile* V = integration -> get_Profile(3,0);

  cout << "I " << I << endl;
  cout << "V " << V << endl << endl;

  cout << "Get_amps" << endl;

  float* I_amps = I -> get_amps();
  float* V_amps = V -> get_amps();

  cout << "Begin LCP standard";

  Reference::To<StandardFlux> lcp = new StandardFlux;

  lcp -> set_fit_shift(true);

  Reference::To<Profile> lcp_profile = new Profile(nbin);

  float* lcp_amps = lcp_profile -> get_amps();

  for(unsigned i = 0; i < nbin ; i++){
    lcp_amps[i] = 0.5*(I_amps [i] -  V_amps[i]);
    cout << lcp_amps[i] << " ";
  }
  cout << "\n\nLoop Finished. Set standard lcp_profile and set flux method to ds[0]\n";

  lcp -> set_standard(lcp_profile);
  ds[0].set_flux_method(lcp);

  cout << "Begin RCP standard\n";

  Reference::To<StandardFlux> rcp = new StandardFlux;

  rcp -> set_fit_shift(true);

  Reference::To<Profile> rcp_profile = new Profile(nbin);

  float* rcp_amps = rcp_profile -> get_amps();

  for(unsigned i = 0; i < nbin; i++){
    rcp_amps[i]  = 0.5*(I_amps [i] + V_amps[i]);
    cout << rcp_amps[i] << " ";
  }

  rcp -> set_standard(rcp_profile);
  ds[1].set_flux_method(rcp);

  }else{
  // Convert
  stdarch = arch->total();
  stdarch->convert_state(Signal::Intensity);

  // Set up DS calculation
  Reference::To<StandardFlux> flux = new StandardFlux;
  flux->set_fit_shift(true); // always init with true, choose later
  flux->set_standard(stdarch->get_Profile(0,0,0));
  ds[0].set_flux_method(flux);
  }
}

void psrflux::setup ()
{

  // If a standard was given, load it.  Otherwise compute a standard
  // from each input file (with warning).
  if (stdfile != "") {
    Reference::To<Archive> arch = Archive::load(stdfile);
    set_standard(arch);
  } else {
    cerr << "Error psrflux: No standard given, will \"self-standard\" each file." << endl;
  }

}

void psrflux::process (Pulsar::Archive* archive)
{

  if(separate_cp_spectra){
    cout << "Computing LCP and RCP spectra. Process method.\n";
    set_standard(archive);

    //compute
    bool single_profile = archive -> get_nsubint()== 1 && archive->get_nchan()==1;
    StandardFlux *flux = dynamic_cast<StandardFlux*>(ds[0].get_flux_method().get());

    //fit shifts
    if(noalign){
      cout<<"No align";
    }else if(align==false && single_profile==false ){
      cout<<"Elif\n";
    }else{
      cout << "set_fit_shift";
      flux->set_fit_shift(true);
    }
    //

    cout << "Finsih";
    
    archive -> convert_state(Signal::Stokes);

    const unsigned nbin = archive -> get_nbin();
    cout << "nbin" << nbin << endl;

    Reference::To<Archive> lcp_archive = archive -> clone();

    lcp_archive -> pscrunch();

    Reference::To<Archive> rcp_archive = archive -> clone();
    rcp_archive -> pscrunch();

    unsigned nsubint = archive -> get_nsubint();
    unsigned nchan = archive -> get_nchan();
   

    //for every channel of every subint
    for(unsigned isub = 0; isub < nsubint; isub++){
      for(unsigned ichan = 0; ichan < nchan; ichan++){
        cout << "isub: " << isub << " chan:" <<ichan << endl;

        //declare pointers
        Profile* lcp_profile = lcp_archive->get_Profile(isub, 0, ichan);
        Profile* rcp_profile = rcp_archive->get_Profile(isub, 0, ichan);

        Profile* I_profile = archive -> get_Profile(isub, 0, ichan);
        Profile* V_profile = archive -> get_Profile(isub, 3, ichan);

        float* I_amps = I_profile -> get_amps();
        float* V_amps = V_profile -> get_amps();

        float* lcp_amps = lcp_profile -> get_amps();
	float* rcp_amps = rcp_profile -> get_amps();

	for(unsigned i = 0; i < nbin; i++){
	lcp_amps[i] = 0.5*(I_amps[i] - V_amps[i]);
	rcp_amps[i] = 0.5*(I_amps[i] + V_amps[i]);
  	//cout << "lcp_amps: " << lcp_amps[i] << " ";
 	//cout << "rcp_amps: " << rcp_amps[i] << endl;
	}

      }
    }

  cout << "Setting Archive: " << endl;

  ds[0].set_Archive(archive);
  ds[1].set_Archive(archive);

  cout << "Compute:" << endl;

  ds[0].compute();
  ds[1].compute();

  cout << "Saving File. Writibg out:" << endl;
  std::string outf = "test_output_2/lcp" + archive->get_filename() + "." + ext;
  cerr << "psrflux: unloading " << outf << endl;

  ds[0].unload(outf, command);

  outf = "test_output_2/rcp" + archive->get_filename() + "." + ext;
  cerr << "psrflux: unloading " << outf << endl;

  ds[1].unload(outf, command);  

  }else{
  // Convert to total intensity
  archive->convert_state(Signal::Intensity);

  // Set self-standard if needed
  if (stdfile=="") set_standard(archive);

  // Test for single-profile data
  bool single_profile = archive->get_nsubint()==1 && archive->get_nchan()==1;

  // Access to the flux computation
  StandardFlux *flux = dynamic_cast<StandardFlux*>(ds[0].get_flux_method().get());

  // If shifts not fit, need to dedisperse and possibly align total
  // with standard.
  if (noalign) {
    archive->dedisperse();
    flux->set_fit_shift(false);
  } else if (align==false && single_profile==false) {
    archive->dedisperse();
    Reference::To<Archive> arch_tot = archive->total();
    Estimate<double> shift = 
      arch_tot->get_Profile(0,0,0)->shift(stdarch->get_Profile(0,0,0));
    stdarch->get_Profile(0,0,0)->rotate_phase(-1.0*shift.get_value());
    flux->set_fit_shift(false);
  } else {
    flux->set_fit_shift(true);
  }

  // Compute DS
  ds[0].set_Archive(archive);
  ds[0].compute();

  // Unload archive with .sm extension
  std::string outf = archive->get_filename() + "." + ext;
  cerr << "psrflux: unloading " << outf << endl;
  ds[0].unload(outf, command);
  }
}

static psrflux program;

int main (int argc, char** argv)
{
  return program.main (argc, argv);
}

