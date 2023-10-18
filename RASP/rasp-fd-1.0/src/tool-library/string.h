#include <string>
#include <vector>

namespace tnstring {

	char const blankchars[] = " \t\n\r";

	//Returns a string with leading/trailing characters of a set stripped
	std::string trim(std::string const&, char const* sepset = blankchars);
	//The same as trim, but returning a right-trimmed string
	std::string rtrim(std::string const& str, char const* sepset = blankchars);
	//The same as trim, but returning a left-trimmed string
	std::string ltrim(std::string const& str, char const* sepset = blankchars);
	//Explode
	std::vector<std::string> explode(const std::string &,const std::string &);
	//Conversion: string to int
	int string2int(const std::string &);
	//Conversion: string to float
	float string2float(const std::string &);
	//Conversion: int to string
	std::string int2string(int);
	//Conversion: float to string
	std::string float2string(float);
	//String comparison, not case sensitive
	bool compare_nocase(std::string, std::string);
	//Conversion: tolower
	std::string toLower(std::string);
};

