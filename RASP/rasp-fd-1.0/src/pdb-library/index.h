//#include <cassert>
#include <iostream>

namespace tnpdb {

	template <class T> class index {
		public:
			// Constructors
			index(int i): _i(i){} 			         // Construct an index from a non-zero integer	
			index():_i(-999999){}                                 // Construct an invalid index
    		// Conversion functions
			operator int() const;                            // Conversion function, for use e.g.: int a = int(index obj), int b = (int)index obj, int c = d + obj
			operator bool() const;                           // Conversion function, for use e.g.: if(obj), if(bool(obj)), if((bool)obj)
			// Overloadiing operator functions
			bool operator == (const index<T> &) const;
			bool operator != (const index<T> &) const;
			bool operator < (const index<T> &) const;
			bool operator > (const index<T> &) const;
			bool operator <= (const index<T> &) const;
			bool operator >= (const index<T> &) const;
			// Methods
			std::ostream &write(std::ostream &) const;
		protected:
			// Members
    		int _i;
	};

	// Definitions
	
	template <class T> index<T>::operator int() const {
		//assert(*this);
		return _i;
	}

	template <class T> index<T>::operator bool() const {
		return (_i != -999999);
	}

	template <class T> bool index<T>::operator == (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i == i._i);
		}
	}
	
	template <class T> bool index<T>::operator != (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i != i._i);
		}
	}

	template <class T> bool index<T>::operator < (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i < i._i);
		}
	}

	template <class T> bool index<T>::operator > (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i > i._i);
		}
	}
	
	template <class T> bool index<T>::operator <= (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i <= i._i);
		}
	}

	template <class T> bool index<T>::operator >= (const index<T> &i) const {
		if (!i || !*this) {
			return false;
		} else {
			return (_i >= i._i);
		}
	}

	template <class T> std::ostream &index<T>::write(std::ostream &out) const {
		if (operator bool()) {
			out << _i;
		} else {
			out << "null";
		}
		return out;
	}

	template <class T> std::ostream &operator << (std::ostream &out, const index<T> &i) {
		return i.write(out);
	}


};

