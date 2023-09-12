#ifndef TRAIN_H
#define TRAIN_H

#include "DistAnal.h"
#include "DihAnal.h"

class Train {
 public:
  Train();
  Train(char *);
  void dist(char *);
  void dih(char *);
  
 private:
  int cutoff;
  double bin;
  string par_dist;
  string par_dih;
};



#endif








