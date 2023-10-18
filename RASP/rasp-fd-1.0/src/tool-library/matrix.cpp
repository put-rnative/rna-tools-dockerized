#include <iostream>
#include <fstream>
#include <cctype>
#include <cstdarg>
#include "matrix.h"

using namespace std;

matrix::~matrix(){
}

matrix::matrix(int _dim,...){
	va_list 
		args;
	int
		_lim,
		len = 1,
		i;
	
	va_start(args,_dim);
	
	dim = _dim;

	//F[0] = 1
	ftr.push_back(1);

	_lim = va_arg(args,int);

	for(i = 0; i < dim && _lim; i++){
		lim.push_back(_lim);
		//F[1] = F[0]*lim[0]; F[2] = F[1]*lim[1]; F[n] = F[n-1]*lim[n-1]
		if(i > 0) ftr.push_back(lim[i-1]*ftr[i-1]);
		_lim = va_arg(args,int);
	}
	va_end(args);
}

void matrix::put(int x,...){
	va_list
		args;

	double
		val;
	int
		map,
		err = 0,
		y,i;
	
	va_start(args,x);
	
	map = 0;// map = F[0]*x + F[1]*y + F[2]*z + ...

	map = map + x*ftr[0];
	if(x >= lim[0]) err = 1;

	for(i = 1; i < dim; i++){
		y = va_arg(args,int);
		map = map + y*ftr[i];
		if(y >= lim[i]) err = 1;
	}
	val = va_arg(args,double);
	va_end(args);

	if(err){
		cerr << "Error: Matrix overflow; no data stored" << endl;
		exit(1);
	}

	if(mtx.size() <= map){
		for(i = mtx.size(); i <= map; i++) mtx.push_back(DEFAULT_V);
	}
	mtx[map] = val;

}

double matrix::get(int x,...){
	va_list
		args;

	int
		map,
		err = 0,
		y,i;

	va_start(args,x);
	
	map = 0;// map = F[0]*x + F[1]*y + F[2]*z + ...
	
	map = map + x*ftr[0];
	if(x >= lim[0]) err = 1;

	for(i = 1; i < dim; i++){
		y = va_arg(args,int);
		map = map + y*ftr[i];
		if(y >= lim[i]) err = 1;
	}
	va_end(args);

	if(err){
		cerr << "Error: Matrix overflow; no data gotten" << endl;
		exit(1);
	}
	
	if(mtx.size() <= map) return DEFAULT_V;
	return mtx[map];
}

const std::vector<int> matrix::get_dimensions() const {
	return lim;
}
