/**********************************************************************
  findSpaceGroup.cpp - Using the spglib to determine the spg of a crystal.

  Copyright (C) 2015 - 2016 by Patrick S. Avery

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 ***********************************************************************/

#include <iostream>

#include "crystal.h"

extern "C" {
#include "../spglib/spglib.h"
}

using namespace std;

int findSpaceGroup(Crystal c, double prec = 0.05);

int main(int argc, char* argv[])
{
  if (argc != 2 && argc != 3) {
    cout << "Usage: ./findSpg <POSCARFileName>\n";
    return -1;
  }

  Crystal c(argv[1]);
  c.wrapAtomsToCell();
  c.printCrystalInfo();

  double tol = 0.05;
  if (argc == 3) tol = atof(argv[2]);

  cout << "Using a tolerance of " << tol << "\n";

  int spg = findSpaceGroup(c);

  cout << "spg is " << spg << "\n";
  return 0;
}

int findSpaceGroup(Crystal c, double prec) {
  // Check that the precision is reasonable
  if (prec < 1e-5) {
    std::cerr  << "findSpaceGroup called with a precision of "
               << prec << ". This is likely an error. Resetting prec to "
               << 0.05 << ".";
    prec = 0.05;
  }

  uint m_spgNumber = 0;
  string m_spgSymbol = "Unknown";
  int num = c.numAtoms();

  vector<vector<double>> latticeVecs = c.getLatticeVecs();

  // if no unit cell or atoms, exit
  if (latticeVecs.size() == 0 || num == 0) {
    std::cerr << "findSpaceGroup( " << prec << " ) called on atom with no cell or atoms!";
    return 0;
  }

  double lattice[3][3];
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      lattice[j][i] = latticeVecs.at(i).at(j);
    }
  }

  // Get atom info
  double (*positions)[3] = new double[num][3];
  int* types = new int[num];

  vector<atomStruct> atoms = c.getAtoms();

  for (int i = 0; i < atoms.size(); i++) {
    types[i]          = atoms.at(i).atomicNum;
    positions[i][0]   = atoms.at(i).x;
    positions[i][1]   = atoms.at(i).y;
    positions[i][2]   = atoms.at(i).z;
  }

  // find spacegroup
  char symbol[21];
  int spg = spg_get_international(symbol,
                                  lattice,
                                  positions,
                                  types,
                                  num, prec);

  delete [] positions;
  delete [] types;

  return spg;
}
