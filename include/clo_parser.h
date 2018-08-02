// Class for parsing command line arguments

#ifndef CLO_PARSER_H
#define CLO_PARSER_H

#include <string>
#include <vector>

/*******************************************************************************

Class to read and process command-line arguments

*******************************************************************************/
class CLOParser {

	private:

	std::vector<std::string> _argv_str;
	std::string _input_file;
	
	/* Converts CLOs to vector of strings */
	
	void readCLOs ( int argc, char *argv[] );

	public:

	/* Constructor */
	
	CLOParser();
	
	/* Checks CLOs for errors */
	
	int checkCLOs (int argc, char *argv[], const std::string & version );
	
	/* Prints various information messages */
	
	void printVersion ( const std::string & version ) const;
	void printHelp () const;
	
	/* Query inputs */
	
	bool requestInputFile () const;
	const std::string & inputFile () const;
};

#endif
