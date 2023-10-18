#ifndef R2D_H
#define R2D_H

#include "../Tools.h"

namespace ss {

struct Nucleotide {
	Nucleotide () {
		seq_ = 'X';
		ss_ = '.';
		num_ = -1;
	}
	Nucleotide (char seq, char ss, int num) {
		seq_ = seq;
		ss_ = ss;
		num_ = num;
	}
	char seq_;
	char ss_;
	int num_;
};

struct RNA_Module {
	RNA_Module () {
		type_ = -1;
		len_ = 0;
		seq_ = "";
		ss_ = "";
		num_ = NULL;
		flag_ = -1;
		son_ = NULL;
		brother_ = NULL;
	}
	RNA_Module (int, int, string, string, int *, int = -1);
	~RNA_Module();
	int type_;
	int len_;
	string seq_;
	string ss_;
	int *num_;
	int flag_;
	RNA_Module *son_;
	RNA_Module *brother_;
};

class R2D {
public:
	R2D(string, string, int = 0);
	~R2D();

	static RNA_Module *setTree(string, string);
	void print();

	RNA_Module *head_;
	RNA_Module *pseudo_head_;

private:
	void delTree(RNA_Module *);
	void printTree(RNA_Module *);
	int len_;
	string seq_;
	string ss_;
	int view_;
};

}; // namespace

#endif






