#include <iostream>
#include <cmath>

#define PI acos(-1.0)

namespace tnpdb {

	class point {
		public:
			inline point(): _x(0), _y(0), _z(0){}                                 //Default
			inline point(double x, double y, double z): _x(x), _y(y), _z(z){}     //Construct a point
			inline double x() const;
			inline double y() const;
			inline double z() const;
			inline double operator [] (unsigned int) const;
			inline point operator = (const point &);
			inline point operator + (const point &) const;
			inline point operator - (const point &) const;
			inline point operator / (double) const;
			inline point operator ^ (const point &) const;
			inline double operator * (const point &) const;
			inline bool operator == (const point &) const;
			inline bool operator != (const point &) const;
		protected:
			double
				_x, _y, _z;
	};

	inline std::ostream &operator << (std::ostream &, const point &); 
	inline std::istream &operator >> (std::istream &, point &);
	inline double distance(const point &, const point &);

	//FUNCTION DEFINITINIOS

	// x
	inline double point::x() const {
		return _x;
	}

	// y
	inline double point::y() const {
		return _y;
	}

	// z
	inline double point::z() const {
		return _z;
	}

	// Get a coordinate by index
	inline double point::operator [] (unsigned int i) const {
		switch(i) {
			case 0: return _x;
			case 1: return _y;
			case 2: return _z;
			default:
				//assert(i < 3);
				return -1;
		}
    }

	// Assignment
    inline point point::operator = (const point &p) {
		_x = p._x;
		_y = p._y;
		_z = p._z;
		return *this;
	}

	// Get the sum between two points
	inline point point::operator + (const point &p) const {
		return point(_x + p._x, _y + p._y, _z + p._z);
	}

	// Get the difference of two points
	inline point point::operator - (const point &p) const {
		return point(_x - p._x, _y - p._y, _z + p._z);
	}

	//Get  division by a scalar
	inline point point::operator / (double scalar) const {
		return point(_x/scalar, _y/scalar, _z/scalar);
	}

	// Get the dot product
    inline double point::operator * (const point &p) const {
		return (_x*p._x + _y*p._y + _z*p._z);
	}

	// Get the cross product
	inline point point::operator ^ (const point &p) const {
		return point(_y*p._z - _z*p._y, _z*p._x - _x*p._z, _x*p._y - _y*p._x);
	}

	// Put a point into a stream
	inline std::ostream &operator << (std::ostream &out, const point &p) {
		out << p.x() << " " << p.y() << " " << p.z();
		return out;
	}

	// Read a point from a stream
	inline std::istream &operator >> (std::istream &in, point &p) {
		double
			x, y, z;
		
		if (in){
			in >> x >> y >> z;
    		p = point(x, y, z);
		}
		
		return in;
	}

	// Check equality
	inline bool point::operator == (const point &p) const {
		if (_x == p._x && _y == p._y && _z == p._z) return true;
		else return false;
	}

	// Check not-equality
	inline bool point::operator != (const point &p) const {
		if (_x != p._x || _y != p._y || _z != p._z) return true;
		else return false;
	}


	// Get the distance between two points
	inline double distance(const point &a, const point &b) {
		return sqrt(((a.x() - b.x())*(a.x() - b.x())) + ((a.y() - b.y())*(a.y() - b.y())) + ((a.z() - b.z())*(a.z() - b.z())));
	}

	// Get the angle between three points
	inline double angle(const point &a, const point &b, const point &c) {
		return acos(((distance(a,b)*distance(a,b)) + (distance(b,c)*distance(b,c)) - (distance(a,c)*distance(a,c)))/(2.0*distance(a,b)*distance(b,c)));
	}
	// Get the arg between three points
	inline double arg(const point &a, const point &b, const point &c) {
		return (((distance(a,b)*distance(a,b)) + (distance(b,c)*distance(b,c)) - (distance(a,c)*distance(a,c)))/(2.0*distance(a,b)*distance(b,c)));
	}
	//Transformation matrix
	inline point transform(const point &a, const double * b[]) {
		return point(a.x()*b[0][0] + a.y()*b[0][1] + a.z()*b[0][2], a.x()*b[1][0] + a.y()*b[1][1] + a.z()*b[1][2], a.x()*b[2][0] + a.y()*b[2][1] + a.z()*b[2][2]);
	}

};

