#include <vector>
#include <string>
#include <cstdlib>

#define DEFAULT_V 0.0

class matrix {
	public:
		matrix(){};
		matrix(int,...);            // Constructor
		~matrix();                  // Destructor
		void put(int,...);
		double get(int,...);
		const std::vector<int> get_dimensions() const;
	private:
		int dim;
		std::vector<int> lim;
		std::vector<int> ftr;
		std::vector<double> mtx;    // Matrix
		
};



