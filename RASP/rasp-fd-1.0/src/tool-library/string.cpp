#include <iostream>
#include <sstream>
#include "string.h"

namespace tnstring {

	//Returns a string with leading/trailing characters of a set stripped
	std::string trim(std::string const& str, char const* sepset){
		std::string::size_type const first = str.find_first_not_of(sepset);
		if(first == std::string::npos){
			return std::string();
		}else{
			return str.substr(first, str.find_last_not_of(sepset) - first + 1);
		}
	}

	//The same as trim, but returning a right-trimmed string
	std::string rtrim(std::string const& str, char const* sepset){
		std::string::size_type const last = str.find_last_not_of(sepset);
		if(last == std::string::npos){
			return std::string();
		}else{
			return str.substr(0, last + 1);
		}
	}
	
	//The same as trim, but returning a left-trimmed string
	std::string ltrim(std::string const& str, char const* sepset){
		std::string::size_type const first = str.find_first_not_of(sepset);
		if(first == std::string::npos){
			return std::string();
		}else{
			return str.substr(first);
		}
	}

	//Explode
	std::vector<std::string> explode(const std::string &sep, const std::string &input) {
		std::vector<std::string> out;
		std::string::size_type
			substring_start(0),
			separator_index(0);
		
		while(std::string::npos != separator_index) {
			separator_index = input.find(sep, substring_start);
			out.push_back(input.substr(substring_start, separator_index - substring_start));
			substring_start = separator_index + sep.length();
		}
		return out;
	}

	//Conversion: string to int
	int string2int(const std::string &str) {
		std::stringstream ss(str);
		int n;
		ss >> n;
		return n;
	}

	//Conversion: string to float
	float string2float(const std::string &str) {
		std::stringstream ss(str);
		float n;
		ss >> n;
		return n;
	}

	//Conversion: int to string
	std::string int2string(int i) {
		std::stringstream ss;
		ss << i;
		return ss.str();
	}
	
	//Conversion: float to string
	std::string float2string(float i) {
		std::stringstream ss;
		ss << i;
		return ss.str();
	}

	//String comparison, not case sensitive
	bool compare_nocase(std::string first, std::string second) {
		unsigned int i = 0;
		while (i < first.length() && i < second.length()) {
			if (tolower(first[i]) < tolower(second[i])) return true;
			else if (tolower(first[i]) > tolower(second[i])) return false;
			++i;
		}
		if (first.length() < second.length()) return true;
		else return false;
	}

	//Conversion: tolower
	std::string toLower(std::string s) {
		//std::stringstream ss;
		for (unsigned int i = 0; i < s.length(); ++i)
		    s[i] = tolower(s[i]);
		return s;
	}

};

