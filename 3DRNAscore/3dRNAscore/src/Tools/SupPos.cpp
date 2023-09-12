#include "SupPos.h"

SupPos::SupPos(Matr_ *m, Matr_ *n) {
	if (m->row != n->row || m->col != 3 || n->col != 3) {
		cerr << "SupPos::SupPos error!" << endl;
		exit(1);
	}
	len = m->row;
	MatrixXf x(len, 3);
	MatrixXf y(len, 3);
	for (int i = 0; i < 3; i++) {
		c1[i] = 0;
		c2[i] = 0;
	}
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 3; j++) {
			c1[j] += m->data[i][j];
			c2[j] += n->data[i][j];
		}
	}
	for (int i = 0; i < 3; i++) {
		c1[i] = c1[i] / len;
		c2[i] = c2[i] / len;
	}
	for (int i = 0; i < len; i++) {
		for (int j = 0; j < 3; j++) {
			x(i, j) = m->data[i][j] - c1[j];
			y(i, j) = n->data[i][j] - c2[j];
		}
	}
	Matrix3f g = x.transpose() * y;

	JacobiSVD<Matrix3f> svd(g, ComputeFullU|ComputeFullV);
	Matrix3f u = svd.matrixU();
	Matrix3f v = svd.matrixV();

	double det = g.determinant();
	Matrix3f I;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) {
				I(i, j) = 1;
			} else {
				I(i, j) = 0;
			}
		}
	}
	if (det < 0) {
		I(2, 2) = -1;
	}
	Matrix3f rotate = u * I * v.transpose();
	
	MatrixXf x_(len, 3);
	x_ = x * rotate;
	x_ = x_ - y;
	rmsd = 0;
	for (int i = 0; i < len; i++) {
		rmsd += x_(i, 0) * x_(i, 0) + x_(i, 1) * x_(i, 1) + x_(i, 2) * x_(i, 2);
	}
	rmsd /= (double) len;
	rmsd = sqrt(rmsd);

	rot = new Matr_(3, 3);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			rot->data[i][j] = rotate(i, j);
		}
	}
}

Matr_ *SupPos::getRot() {
	Matr_ *r = new Matr_(rot);
	return r;
}

Point *SupPos::getC1() {
	Point *p = new Point(c1);
	return p;
}

Point *SupPos::getC2() {
	Point *p = new Point(c2);
	return p;
}

double SupPos::getRMSD() {
	return rmsd;
}

SupPos::~SupPos() {
	delete rot;
}








