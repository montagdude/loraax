#include <string>
#include <vector>
#include <iostream>
#include "util.h"
#include "clo_parser.h"

/*******************************************************************************

Reads command line arguments and saves as vector of strings

*******************************************************************************/
void CLOParser::readCLOs ( int argc, char *argv[] )
{
	int i;
	
	_argv_str.resize(0);
	for ( i = 0; i < argc; i++ )
	{
		_argv_str.push_back(chararray2string(argv[i]));
	}
}

/*******************************************************************************

Constructor

*******************************************************************************/
CLOParser::CLOParser ()
{
	_argv_str.resize(0);
	_input_file = "";
}

/*******************************************************************************

Checks CLOs for errors. Returns 0 on successful read, 1 if there are errors, or
-1 if the user asked to print help or version info.

*******************************************************************************/
int CLOParser::checkCLOs ( int argc, char *argv[], const std::string & version )
{
	int i;
	
	// Save CLOs as vector of strings
	
	readCLOs(argc, argv);
	
	// Check for errors
	
	if (argc < 2)
	{
		printHelp();
		return 1;
	}

	// Check for --help or --version

	i = 1;
	while (i < argc)
	{
		if ( (_argv_str[i] == "--help") || (_argv_str[i] == "-h") )
		{
			printHelp();
			return -1;
		}
		else if ( (_argv_str[i] == "--version") || (_argv_str[i] == "-v") )
		{
			printVersion(version);
			return -1;
		}
		i += 1;
	}

	// Read input file

	_input_file = _argv_str[1];

	return 0;
}

/*******************************************************************************

Prints various information messages

*******************************************************************************/
void CLOParser::printVersion ( const std::string & version ) const
{
	std::cout << "LOw Reynolds number Aircraft Analysis with Xfoil, version "
	          << version << std::endl;
	std::cout << "Copyright (C) 2019 Daniel Prosser" << std::endl;
	std::cout << "License GPLv2: GNU GPL version 2:\n"
	          << "<https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>"
	          << std::endl;
	std::cout << "This is free software; you are free to change it and "
	          << "redistribute it." << std::endl;
	std::cout << "There is NO WARRANTY, to the extent permitted by law."
	          << std::endl;
}

void CLOParser::printHelp () const
{
	std::cout << "Usage: loraax [-h] [-v] INPUTFILE" << std::endl;
	std::cout << std::endl;
	std::cout << "Required inputs:" << std::endl;
	std::cout << "INPUTFILE                Perform analysis defined in XML "
	          << "input file" << std::endl;
	std::cout << "  -h, --help             Display usage information and exit"
	          << std::endl;
	std::cout << "  -v, --version          Display version number and exit"
	          << std::endl;
	std::cout << std::endl;
	std::cout << "loraax home page: https://github.com/montagdude/loraax"
	          << std::endl;
}

/*******************************************************************************

Query inputs

*******************************************************************************/
bool CLOParser::requestInputFile () const
{
	if (_input_file == "") { return false; }
	else { return true; }
}

const std::string & CLOParser::inputFile () const { return _input_file; }
