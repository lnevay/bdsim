#include "BDSExecOptions.hh"

#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <string>

#include "globals.hh" // geant4 types / globals

#include "BDSDebug.hh"
#include "BDSMaterials.hh"
#include "BDSOutputLoader.hh"
#include "BDSUtilities.hh"

#include "parser/getEnv.h"
#include "parser/options.h"

BDSExecOptions::BDSExecOptions(int argc, char **argv):
  options(GMAD::Options()),
  ignoreSIGINT(false)
{
  Parse(argc, argv);
  /// after parsing the absolute path can be reconstructed  
  options.set_value("bdsimPath", std::string(GetPath(options.inputFileName)));
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "BDSIMPATH set to: " << options.bdsimPath << G4endl;
#endif
  if (options.recreate)
    {
      BDSOutputLoader loader(options.recreateFileName);
      GMAD::Options recreateOptions = loader.Options();
      // Give precedence to exec options - only ones that have been set.
      recreateOptions.Amalgamate(options, true);
      options = recreateOptions; // Now replace member.
    }
}

void BDSExecOptions::Parse(int argc, char **argv)
{
  static struct option LongOptions[] = {{ "help" , 0, 0, 0 },
					{ "verbose", 0, 0, 0 },
					{ "verbose_step", 0, 0, 0 },
					{ "verbose_event", 0, 0, 0 },
					{ "verbose_event_num", 1, 0, 0 },
					{ "verbose_G4run", 1, 0, 0 },
					{ "verbose_G4event", 1, 0, 0 },
					{ "verbose_G4tracking", 1, 0, 0 },
					{ "verbose_G4stepping", 1, 0, 0 },
					{ "file", 1, 0, 0 },
					{ "distrFile", 1, 0, 0 },
					{ "vis_debug", 0, 0, 0 },
					{ "vis_mac", 1, 0, 0 },
					{ "output", 1, 0, 0 },
					{ "outfile", 1, 0, 0 },
					{ "batch", 0, 0, 0 },
					{ "materials", 0, 0, 0 },
					{ "circular", 0, 0, 0 },
					{ "seed",           1, 0, 0},
					{ "recreate",       1, 0, 0},
					{ "startFromEvent", 1, 0, 0},
					{ "writeseedstate", 0, 0, 0},
					{ "seedstate",      1, 0, 0},
					{ "seedStateFileName", 1, 0, 0},
					{ "survey", 1, 0, 0 },
					{ "ngenerate", 1, 0, 0 },
					{ "nGenerate", 1, 0, 0 },
					{ "exportgeometryto", 1, 0, 0 },
					{ "generatePrimariesOnly", 0, 0, 0 },
					{ "ignoresigint", 0, 0, 0},
					{ 0, 0, 0, 0 }};
  
  int OptionIndex  = 0;
  int OptionNumber = 0;
  const char* optionName;
  // number conversion check
  bool conversion = true;
  for(;;) {
    OptionIndex = 0;
  
    // see e.g. http://linux.die.net/man/3/getopt
    int c = getopt_long(argc, argv, "Vv",
			LongOptions, &OptionIndex );
    
    if ( c == -1 ) // end of options list
      break;

    switch (c) {
    case '?': // unrecognised option
      G4cout << "invalid option for command " << argv[0] << G4endl << G4endl << G4endl;
      Usage();
      exit(1);
      break;
    case 0:
      OptionNumber++;
      optionName = LongOptions[OptionIndex].name;
      if( !strcmp(optionName , "help") )
	{
	  Usage();
	  exit(0);
	}      
      else if( !strcmp(optionName , "batch") )
	{options.set_value("batch",true);}
      else if( !strcmp(optionName , "verbose") )
	{options.set_value("verbose", true);}
      else if( !strcmp(optionName , "verbose_step") )
	{// we shouldn't have verbose steps without verbose events
	  options.set_value("verboseStep", true); 
	  options.set_value("verboseEvent", true);
	}
      else if( !strcmp(optionName , "verbose_event") )
	{options.set_value("verboseEvent", true);}
      else if( !strcmp(optionName , "verbose_event_num") )
	{
	  int result = -1;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("verboseEventNumber", result);
	}
      else if( !strcmp(optionName , "verbose_G4run") )
	{
	  int result = 0;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("verboseRunLevel", result);
	}
      else if( !strcmp(optionName , "verbose_G4event") )
	{
	  int result = 0;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("verboseEventLevel", result);
	}
      else if( !strcmp(optionName , "verbose_G4tracking") )
	{
	  int result = 0;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("verboseTrackingLevel", result);
	}
      else if( !strcmp(optionName , "verbose_G4stepping") )
	{
	  int result = 0;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("verboseSteppingLevel", result);
	}
      else if( !strcmp(optionName , "output") )
	{options.set_value("outputFormat", std::string(optarg));}
      else if( !strcmp(optionName , "outfile") )
	{options.set_value("outputFileName", std::string(optarg));}
      else if( !strcmp(optionName , "survey") )
	{
	  options.set_value("surveyFileName", std::string(optarg));
	  options.set_value("survey",         true);
	}
      else if( !strcmp(optionName , "file") )
	{options.set_value("inputFileName", std::string(optarg));}
      else if( !strcmp(optionName, "distrFile") )
	{// build absolute path
          options.set_value("distrFile", BDS::GetCurrentDir() + "/" + std::string(optarg));
        }
      else if( !strcmp(optionName , "vis_debug") )
	{options.set_value("visDebug", true);}
      else if( !strcmp(optionName , "vis_mac") )
	{options.set_value("visMacroFileName", BDS::GetCurrentDir() + "/" + std::string(optarg));}
      else if( !strcmp(optionName, "materials") )
	{
	  BDSMaterials::Instance()->ListMaterials();
	  exit(0);// return after printing material list
	}
      else if( !strcmp(optionName, "circular") )
	{options.set_value("circular", true);}
      else if( !strcmp(optionName, "seed")  )
	{
	  int result = -1;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("seed", result);
	}
      else if( !strcmp(optionName, "recreate") )
	{
          options.set_value("recreate", true);
          options.set_value("recreateFileName", BDS::GetCurrentDir() + "/" + std::string(optarg));
	}
      else if( !strcmp(optionName, "startFromEvent") )
	{
	  int result = 0;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("startFromEvent", result);
	}
      else if( !strcmp(optionName, "writeseedstate") )
	{options.set_value("writeSeedState", true);}
      else if( !strcmp(optionName, "seedstate")  || !strcmp(optionName, "seedStateFileName"))
	{
          options.set_value("useASCIISeedState", true);
          options.set_value("seedStateFileName", BDS::GetCurrentDir() + "/" + std::string(optarg));
	}
      else if( !strcmp(optionName, "ngenerate") || !strcmp(optionName, "nGenerate"))
	{
	  int result = 1;
	  conversion = BDS::IsInteger(optarg, result);
	  options.set_value("ngenerate", result);
	  options.set_value("matchDistribFileLength", false); // ngenerate overrides.
	}
      else if( !strcmp(optionName, "generatePrimariesOnly") )
	{options.set_value("generatePrimariesOnly", true);}
      else if( !strcmp(optionName, "ignoresigint") )
	{ignoreSIGINT = true;}
      else if( !strcmp(optionName, "exportgeometryto") )
	{// TBC - this should be put into geometry classes
	  std::string fn = optarg;
	  if (fn.substr(fn.find_last_of(".") + 1) == "gdml")
	    {
	      options.set_value("exportType",     std::string("gdml"));
	      options.set_value("exportFileName", fn);
	    }
	  else
	    {
	      G4cerr << __METHOD_NAME__ << "Unknown geometry format \""
		     << fn.substr(fn.find_last_of(".") + 1) << "\"\n"
		     << "Please specify a valid filename extension - options are: \"gdml\"" << G4endl;
	      exit(1);
	    }
	  options.set_value("exportGeometry", true);
	}
      // remember if you extend this to do it also in the usage print out
	  
      if (conversion == false)
	{
	  // conversion from character string to number went wrong, exit
	  G4cerr << __METHOD_NAME__ << "Conversion to number (or integer) went wrong for \""
		 << optionName << "\" with value: \"" << optarg << "\"" << G4endl;
	  exit(1);
	}

      break;
      
    default:
      G4cout << "Warning unknown returned character code " <<  c << G4endl;
      break;
    }
  }
  // there should be no remaining options
  if (OptionNumber < argc - 1)
    {
      G4cout << "ERROR there are remaining options: " << G4endl;
      for (int i=1; i<argc; i++)
	{
	  // options with '-' are ignored by getopt_long, other unknown options are covered
	  if (strncmp(argv[i], "-", 1))
	    {G4cout << argv[i] << G4endl;}
	}
      
      G4cout << "Please check your input" << G4endl;
      exit(1);
    }
  }

void BDSExecOptions::Usage() const
{
  G4cout<<"Usage: bdsim [options]"<<G4endl;
  G4cout<<"Options:"<<G4endl;
  G4cout<<"--file=<filename>         : specify the lattice and options file "<<G4endl
	<<"--batch                   : batch mode - no graphics"<<G4endl
	<<"--circular                : assume circular machine - turn control"<<G4endl
	<<"--distrFile=<file>        : specify which file to use for the bunch distribution" << G4endl
	<<"--exportgeometryto=<file> : export the geometry to a file - extension determines format"<<G4endl
	<<"                            where possible extensions are (\"gdml\")"<<G4endl
	<<"--generatePrimariesOnly   : generate N primary particle coordinates without simulation then quit"<<G4endl
	<<"--help                    : display this message"<<G4endl
	<<"--materials               : list materials included in bdsim by default"<<G4endl
	<<"--output=<fmt>            : output format (rootevent|ascii|combined|none), default rootevent"<<G4endl
	<<"--outfile=<file>          : output file name. Will be appended with _N"<<G4endl
        <<"                            where N = 0, 1, 2, 3... etc."<<G4endl
	<<"--ngenerate=N             : the number of primary events to simulate - overrides the ngenerate " << G4endl
	<<"                            option in the input gmad file" << G4endl
        <<"--seed=N                  : the seed to use for the random number generator" << G4endl
	<<"--seedStateFileName=<file>: use this ASCII file seed state to run an event" << G4endl
	<<"--recreate=<file>         : the rootevent file to recreate events from" << G4endl
	<<"--startFromEvent=N        : event offset to start from when recreating events" << G4endl
	<<"--survey=<file>           : print survey info to <file>"<<G4endl
	<<"--verbose                 : display general parameters before run"<<G4endl
	<<"--verbose_event           : display information for every event "<<G4endl
	<<"--verbose_event_num=N     : display tracking information for event number N"<<G4endl
	<<"--verbose_step            : display tracking information after each step"<<G4endl
	<<"--verbose_G4event=N       : set Geant4 Event manager verbosity level"<<G4endl
	<<"--verbose_G4run=N         : set Geant4 verbosity level (see Geant4 manual for details)"<<G4endl
	<<"--verbose_G4stepping=N    : set Geant4 Stepping manager verbosity level"<<G4endl
	<<"--verbose_G4tracking=N    : set Geant4 Tracking manager verbosity level [-1:5]"<<G4endl
	<<"--vis_debug               : display all volumes in visualiser"<<G4endl
	<<"--vis_mac=<file>          : file with the visualisation macro script, default provided by BDSIM openGL (OGLSQt))" << G4endl
	<<"--writeseedstate          : write an ASCII file seed state for each event" << G4endl;
}

void BDSExecOptions::Print() const
{
  G4cout << std::boolalpha; // print textual representation of bool
  G4cout << __METHOD_NAME__ << std::setw(23) << " inputFileName: "       << std::setw(15) << options.inputFileName       << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " batch: "               << std::setw(15) << options.batch               << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " circular: "            << std::setw(15) << options.circular            << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " distrFile: "           << std::setw(15) << options.distribFile         << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " exportgeometryto "     << std::setw(15) << options.exportFileName      << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " generatePrimariesOnly "<< std::setw(15) << options.generatePrimariesOnly << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " ngenerate: "           << std::setw(15) << options.nGenerate           << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " outputFileName: "      << std::setw(15) << options.outputFileName      << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " outputFormat: "        << std::setw(15) << options.outputFormat        << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " seed: "                << std::setw(15) << options.seed                << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " seedstate: "           << std::setw(15) << options.seedStateFileName   << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " survey: "              << std::setw(15) << options.survey              << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " surveyFileName: "      << std::setw(15) << options.surveyFileName      << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " verbose: "             << std::setw(15) << options.verbose             << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseEvent: "        << std::setw(15) << options.verboseEvent        << G4endl;  
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseStep: "         << std::setw(15) << options.verboseStep         << G4endl;  
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseRunLevel: "     << std::setw(15) << options.verboseRunLevel     << G4endl;  
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseEventLevel: "   << std::setw(15) << options.verboseEventLevel   << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseTrackingLevel: "<< std::setw(15) << options.verboseTrackingLevel<< G4endl;  
  G4cout << __METHOD_NAME__ << std::setw(23) << " verboseSteppingLevel: "<< std::setw(15) << options.verboseSteppingLevel<< G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " visMacroFileName: "    << std::setw(15) << options.visMacroFileName    << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " visDebug: "            << std::setw(15) << options.visDebug            << G4endl;
  G4cout << __METHOD_NAME__ << std::setw(23) << " ignoreSIGINT: "        << std::setw(15) << ignoreSIGINT                << G4endl;
  G4cout << std::noboolalpha; // reset to default printing
}

G4String BDSExecOptions::GetPath(G4String fileName)
{
  // Set fullPath to mirror what is done in parser.l (i.e. if no environment
  // variable set, assume base filename path is that of the gmad file).
  G4String fullPath = getEnv("BDSIMPATH");
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "filename = " << fileName << G4endl;
#endif
  if(fullPath.empty())
    {
      G4String inputFilepath = "";
      // get the path part of the supplied path to the main input file
      G4String::size_type found = fileName.rfind("/"); // find the last '/'
      if (found != G4String::npos)
	{inputFilepath = fileName.substr(0,found);} // the path is the bit before that
      // else remains empty string
      // need to know whether it's an absolute or relative path
      if ((fileName.substr(0,1)) == "/")
	{// the main file has an absolute path
	  fullPath = inputFilepath;
	}
      else
	{fullPath = BDS::GetCurrentDir() + "/" + inputFilepath;}
    }
  // add additional slash just to be safe
  fullPath += "/";
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "fullpath = " << fullPath << G4endl;
#endif
  return fullPath;
}
