/*
 * RunPredict.cpp
 *
 *  Created on: 2020Äê12ÔÂ23ÈÕ
 *      Author: pengx
 */

#include "forcefield/NeighborConnectionGenerator.h"
#include "geometry/localframe.h"
#include "model/BaseDistanceMatrix.h"
#include "model/ProteinRep.h"
#include "pred/BRFoldingTree.h"
#include "pred/BRFoldingTreeBasic.h"
#include "pred/MCRun.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace NSPmodel;
using namespace NSPforcefield;
using namespace NSPpred;
using namespace std;

void checkInputFile(const string &inputFile) {
  ifstream file;
  cout << "check input file: " << endl;
  file.open(inputFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open inputFile: " << inputFile << endl;
    exit(0);
  }
  file.close();

  NSPtools::InputParser input(inputFile);
  string pdbFile = input.getValue("pdb");
  string baseSeq = input.getValue("seq");
  string baseSec = input.getValue("sec");
  string nwcSec = input.getValue("nwc");

  if (pdbFile == "") {
    cout << "invalid inputFile: initial pdb required" << endl;
    exit(0);
  }

  if (baseSeq == "") {
    cout << "invalid inputFile: base sequence required" << endl;
    exit(0);
  }

  if (baseSec == "") {
    cout << "invalid inputFile: Watson-Crick pairing sequence required" << endl;
    exit(0);
  }

  if (nwcSec == "") {
    cout << "invalid inputFile: non-Watson-Crick pairing sequence required"
         << endl;
    exit(0);
  }

  if (baseSec.length() != baseSeq.length()) {
    cout << "the length of Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    exit(0);
  }

  if (nwcSec.length() != baseSeq.length()) {
    cout << "the length of non-Watson-Crick pairing sequence must equal to the "
            "length of base sequence"
         << endl;
    exit(0);
  }

  file.open(pdbFile.c_str(), ios::in);
  if (!file.is_open()) {
    cout << "fail to open pdb file: " << pdbFile << endl;
    exit(0);
  }
}

void printHelp() { cout << "Usage: BRiQ_Energy $INPUT_FILE" << endl; }

int main(int argc, char **argv) {

  // Usage: BRiQ_Predict $INPUT_FILE $OUTPUT_PDB $RANDOM_SEED
  if (argc != 2 || argv[1] == "-h") {
    printHelp();
    exit(0);
  }

  srand(time(NULL));

  string inputFile = string(argv[1]);
  checkInputFile(inputFile);

  EnergyTable *et = new EnergyTable();
  BRFoldingTree *ft = new BRFoldingTree(inputFile, et);
  cerr << inputFile << " " << ft->totalEnergy(1.0, 1.0, 1.0, 0.0, false)
       << endl;

  delete et;
  delete ft;
}
