// Contains basic utility functions and constants

#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <vector>

// Public routines

void conditional_stop ( int, const std::string &, const std::string & );
void print_warning ( const std::string &, const std::string & );
void go_to_line ( std::ifstream & , unsigned int );
int string2int ( const std::string & );
double string2double ( const std::string & );
std::string int2string ( int );
std::string double2string ( const double & );
std::string bool2string ( bool );
bool string2bool ( const std::string & );
bool is_space ( char );
std::string get_name_end ( const std::string & );
std::string get_name_begin ( const std::string & );
std::string bracket_name ( const std::string & );
std::string repeat_string ( const std::string & , unsigned int );
std::vector<std::string> split_string ( const std::string & );
std::vector<std::string> split_string ( const std::string &, char );

#endif
