// Contains basic utility functions

#define _USE_MATH_DEFINES

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

/******************************************************************************/
//
// Conditionally stops code execution
//
/******************************************************************************/
void conditional_stop ( int ecode, const std::string & location, 
                        const std::string & message)
{
	// Exit only if exit code != 0
	
	if (ecode != 0)
	{
		std::cerr << "Error in " << location << ": " << message << std::endl;
		exit(EXIT_FAILURE);
	}
}

/******************************************************************************/
//
// Prints warning message
//
/******************************************************************************/
void print_warning ( const std::string & location, const std::string & message)
{
	std::cout << "Warning in " << location << ": " << message << std::endl;
}

/******************************************************************************/
//
// Goes to a specific line in a file
//
/******************************************************************************/
void go_to_line ( std::ifstream & infile, unsigned int linenum )
{
	unsigned int i;
	std::string line;
	
	// Go to the beginning of the file
	
	infile.clear();   // In case eof has been reached
	infile.seekg(0);
	
	// Skip as many lines as requested
	
	for ( i = 0; i < linenum; i++ ) { getline(infile, line); } 
	
	return;
}

/******************************************************************************/
//
// Converts string to int. Returns 0 on success, 1 on failure.
//
/******************************************************************************/
int string2int ( const std::string & instr, int & outval )
{
	std::stringstream ss(instr);
	ss >> outval;
	
	if (! ss)
		return 1;
	else
		return 0;
}

/******************************************************************************/
//
// Converts string to double. Returns 0 on success, 1 on failure.
//
/******************************************************************************/
int string2double ( const std::string & instr, double & outval )
{
	std::stringstream ss(instr);
	ss >> outval;
	
	if (! ss)
		return 1;
	else
		return 0;
}

/******************************************************************************/
//
// Converts int to string
//
/******************************************************************************/
std::string int2string ( int inval )
{
	std::stringstream ss;
	std::string outstr;
	ss << inval;
	return ss.str();
}

/******************************************************************************/
//
// Converts double to string
//
/******************************************************************************/
std::string double2string ( const double & inval )
{
	std::stringstream ss;
	std::string outstr;
	ss << inval;
	return ss.str();
}

/******************************************************************************/
//
// Converts bool to string
//
/******************************************************************************/
std::string bool2string ( bool inval )
{
	std::string outstr;
	outstr = "false";
	if (inval) { outstr = "true"; }
	return outstr;
}

/*******************************************************************************

Converts a char array to string

*******************************************************************************/
std::string chararray2string ( const char chararray[] )
{
	std::string outstr;
	std::stringstream ss;
	
	ss << chararray;
	ss >> outstr;
	return outstr;
}

/******************************************************************************/
//
// Converts string to bool
//
/******************************************************************************/
bool string2bool ( const std::string & instr )
{
	if ( (instr == "TRUE") || (instr == "True") || (instr == "true") )
	{
		return true;
	}
	else { return false; }
}

/******************************************************************************/
//
// Determines whether a char is a space.  A space can either be ' ' or '\000' 
// (trailing spaces).
//
/******************************************************************************/
bool is_space ( char inchar )
{
	if ( inchar == ' ' || inchar == '\000' ) { return true; }
	else { return false; }
}

/******************************************************************************/
//
// Removes white space and comments from the end of a string ("!" denotes 
// comment)
//
/******************************************************************************/
std::string get_name_end ( const std::string & name )
{
	int i, len, comment;
	std::string newstr;
	
	// Length of string
	
	len = name.length();
	
	// Test to see if it contains "!"
	
	comment = name.find("!");
	if ( comment != -1 ) { len = comment; }
	
	// Remove any trailing white space
	
	for ( i = len - 1 ; i >= 0; i-- )
	{
		if ( is_space(name[i]) ) { comment = i; }
		else { break; }
	}
	
	if ( comment != -1 ) { len = comment; }
	
	// Create the output string
	
	newstr = name.substr(0, len);
	
	return newstr;
}

/******************************************************************************/
//
// Removes white space from the beginning of a string
//
/******************************************************************************/
std::string get_name_begin ( const std::string & name )
{
	int i, len, comment;
	std::string newstr;
	
	// Length of string
	
	len = name.length();
	
	// Remove any leading white space
	
	comment = 0;
	for ( i = 0; i < len; i++ )
	{
		if ( is_space(name[i]) ) { comment = i + 1; }
		else { break; }
	}
	
	// Create the output string
	
	newstr = name.substr(comment, len);
	
	return newstr;
}

/******************************************************************************/
//
// Removes white space and comments from the end of a string, and then removes 
// white space from beginning of the string.
//
/******************************************************************************/
std::string bracket_name ( const std::string & name )
{
	std::string str1, str2;
	
	// Remove trailing white space and comments from the end
	
	str1 = get_name_end(name);
	
	// Remove white space from the beginning
	
	str2 = get_name_begin(str1);
	
	return str2;
}

/******************************************************************************/
//
// Returns input string copied and appended nrepeat times
//
/******************************************************************************/
std::string repeat_string ( const std::string & instr, unsigned int nrepeat )
{ 
	std::string outstr;
	unsigned int i;
	
	outstr = instr;
	for ( i = 1; i < nrepeat; i++ ) { outstr += instr; }
	return outstr;
}

/******************************************************************************/
//
// Converts string to vector of strings using spaces as the delimiter
//
/******************************************************************************/
std::vector<std::string> split_string ( const std::string & instr )
{
	unsigned int i, len, first, last;
	bool item;
	std::string tempstr;
	std::vector<std::string> outvec;
	
	// Remove leading and trailing white space
	
	tempstr = bracket_name(instr);
	
	// Read items
	
	len = tempstr.length();
	item = false;
	for ( i = 0; i < len; i++ )
	{
		// Previous character was a space
		
		if (not item)
		{
			// New item found
			
			if (not is_space(tempstr[i]))
			{
				item = true;
				first = i;
			}
		}
		
		// Previous character was not a space
		
		else
		{
			// End of item found -> add to output vector
			
			if (is_space(tempstr[i]))
			{
				item = false;
				last = i - 1;
				outvec.push_back(tempstr.substr(first, last-first+1));
			}
		}
		
		// Handling end of string
		
		if ( (not is_space(tempstr[i])) and (i == len - 1) and (item) )
		{
			item = false;
			last = i;
			outvec.push_back(tempstr.substr(first, last-first+1));
		}
	}
	
	return outvec;
}

/******************************************************************************/
//
// Converts string to vector of strings with a specified delimiter
//
/******************************************************************************/
std::vector<std::string> split_string ( const std::string & instr, char delim )
{
	unsigned int i, len, first, last;
	bool item;
	std::string tempstr;
	std::vector<std::string> outvec;
	
	// Remove leading and trailing white space
	
	tempstr = bracket_name(instr);
	
	// Read items
	
	len = tempstr.length();
	item = false;
	for ( i = 0; i < len; i++ )
	{
		// Previous character was a space or delimiter
		
		if (not item)
		{
			// New item found
			
			if (not is_space(tempstr[i]))
			{
				item = true;
				first = i;
			}
		}
		
		// Previous character was not a space
		
		else
		{
			// End of item found -> add to output vector
			
			if (tempstr[i] == delim)
			{
				item = false;
				last = i - 1;
				outvec.push_back(tempstr.substr(first, last-first+1));
			}
		}
		
		// Handling end of string
		
		if ( (tempstr[i] != delim) and (i == len - 1) and (item) )
		{
			item = false;
			last = i;
			outvec.push_back(tempstr.substr(first, last-first+1));
		}
	}
	
	return outvec;
}

/******************************************************************************/
//
// Gets the minimum value of a std::vector<double>
//
/******************************************************************************/
double vector_min ( const std::vector<double> & vec )
{
	unsigned int i, n;
	double fmin;
	
	n = vec.size();
	fmin = 1.E+12;
	for ( i = 0; i < n; i++ )
	{
		if (vec[i] < fmin)
			fmin = vec[i];
	}
	
	return fmin;
}

/******************************************************************************/
//
// Returns -1 for val < 0 and 1 for val >= 0
//
/******************************************************************************/
double sign ( const double & val )
{
	if (val < 0.)
		return -1.;
	else
		return 1.;
}
