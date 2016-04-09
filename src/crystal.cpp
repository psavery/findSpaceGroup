/**********************************************************************
  crystal.cpp - Custom crystal class for generating crystals with
                specific space groups.

  Copyright (C) 2015 - 2016 by Patrick S. Avery

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

 ***********************************************************************/

#include <iostream>
// For setprecision()
#include <iomanip>
// for sqrt(), sin(), cos(), etc.
#include <cmath>
// for writing to POSCAR
#include <fstream>

#include "crystal.h"
#include "utilityFunctions.h"

// For atomic radii and symbols
#include "elemInfo.h"

using namespace std;

//#define CRYSTAL_DEBUG
//#define NEAREST_NEIGHBOR_DEBUG
//#define CENTER_CELL_DEBUG
//#define IAD_DEBUG

Crystal::Crystal(latticeStruct l, vector<atomStruct> a, bool usingVdwRad) :
  m_lattice(l),
  m_atoms(a),
  m_usingVdwRadii(usingVdwRad)
{

}

Crystal::Crystal(const std::string& filename, bool usingVdwRad) :
  m_lattice(latticeStruct()),
  m_atoms(vector<atomStruct>()),
  m_usingVdwRadii(usingVdwRad)
{
  latticeStruct ls;
  vector<atomStruct> as;
  if (!readPOSCAR(filename, ls, as)) {
    cerr << "Failed to construct a crystal after attempting to read the "
         << "POSCAR file: " << filename << "\n";
  }
  else {
    m_lattice = ls;
    m_atoms = as;
  }
}

void Crystal::removeAtomAt(size_t i)
{
  if (i >= m_atoms.size())
    std::cout << "Error: tried to remove an atom at index " << i << " and the "
              << "size is only " << m_atoms.size() << "!\n";
  else m_atoms.erase(m_atoms.begin() + i);
}

void Crystal::removeAtom(const atomStruct& as)
{
  int ind = getAtomIndexNum(as);
  if (ind == -1) {
    cout << "Error: " << __FUNCTION__ << " was called to remove an atom that "
         << "is not a member of the cell!\n";
    return;
  }
  removeAtomAt(ind);
}

void Crystal::removeAllNewAtomsSince(const atomStruct& as)
{
  // Since atoms get appended to the vector in order, we assume all indices
  // including and greater than our current one are new
  for (size_t i = getAtomIndexNum(as) + 1; i < m_atoms.size(); i++) {
    removeAtomAt(i);
    // size will keep going down until we stop
    i--;
  }
}

static inline bool atomsHaveSamePosition(const atomStruct& a1,
                                         const atomStruct& a2)
{
  // min distance
  double minD = 0.00001;
  double dx = fabs(a1.x - a2.x);
  double dy = fabs(a1.y - a2.y);
  double dz = fabs(a1.z - a2.z);
  if (dx < 0.00001 && dy < 0.00001 && dz < 0.00001) return true;
  return false;
}

// We're assuming we're using fractional coordinates...
static inline void wrapUnitToCell(double& u)
{
  double minD = 0.00001;
  while (u < 0.0) u += 1.0;
  while (u >= 1.0 || fabs(u - 1.0) < minD) u -= 1.0;
}

static inline void wrapAtomToCell(atomStruct& as)
{
  wrapUnitToCell(as.x);
  wrapUnitToCell(as.y);
  wrapUnitToCell(as.z);
}

void Crystal::wrapAtomsToCell()
{
  for (size_t i = 0; i < m_atoms.size(); i++) wrapAtomToCell(m_atoms[i]);
}

void Crystal::removeAtomsWithSameCoordinates()
{
  for (size_t i = 0; i < m_atoms.size(); i++) {
    for (size_t j = i + 1; j < m_atoms.size(); j++) {
      if (atomsHaveSamePosition(m_atoms.at(i), m_atoms.at(j))) {
        removeAtomAt(j);
        j--;
      }
    }
  }
}

bool Crystal::addAtomIfPositionIsEmpty(atomStruct& as)
{
  wrapAtomToCell(as);
  bool positionIsEmpty = true;
  for (size_t i = 0; i < m_atoms.size(); i++) {
    if (atomsHaveSamePosition(as, m_atoms.at(i)) &&
        as.atomicNum == m_atoms.at(i).atomicNum) {
      positionIsEmpty = false;
      break;
    }
  }

  if (positionIsEmpty) {
    addAtom(as);
    return true;
  }

  return false;
}

atomStruct Crystal::getAtomInCartCoords(const atomStruct& as) const
{
  atomStruct ret = as;
  const double& a = m_lattice.a;
  const double& b = m_lattice.b;
  const double& c = m_lattice.c;
  const double alpha = deg2rad(m_lattice.alpha);
  const double beta = deg2rad(m_lattice.beta);
  const double gamma = deg2rad(m_lattice.gamma);
  const double v = getUnitVolume();
  // We're just gonna do the matrix multiplication by hand...
  ret.x = a * as.x + b * cos(gamma) * as.y + c * cos(beta) * as.z;
  ret.y = b * sin(gamma) * as.y +
          c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma) * as.z;
  ret.z = c * v / sin(gamma) * as.z;
  // I want to make sure these are all positive
  // When the value is zero, unfortunately, it can be very slightly negative
  // sometimes... So I commented out this assertion
  //assert(ret.x > 0 && ret.y > 0 && ret.z > 0);
  return ret;
}

vector<uint> Crystal::getVectorOfAtomicNums() const
{
  vector<uint> ret;
  for (size_t i = 0; i < m_atoms.size(); i++) {
    ret.push_back(m_atoms.at(i).atomicNum);
  }
  return ret;
}

double Crystal::getUnitVolume() const
{
  const latticeStruct& l = m_lattice;
  const double alpha = deg2rad(m_lattice.alpha);
  const double beta = deg2rad(m_lattice.beta);
  const double gamma = deg2rad(m_lattice.gamma);
  return sqrt(1.0 -
              pow(cos(alpha), 2.0) -
              pow(cos(beta), 2.0) -
              pow(cos(gamma), 2.0) +
              2.0 * cos(alpha) * cos(beta) * cos(gamma));
}

vector<vector<double>> Crystal::getLatticeVecs() const
{
  // To do this, we are going to use a little "hack" using code
  // I've already written
  // For a, we are going to a create a temporary atomStruct object that
  // has fractional coordinates of (1,0,0). Then, we will convert
  // the coordinates to Cartesian using another function I wrote.
  // That will produce an atom whose position is exactly the same as the
  // first vector. Follow this same logic for the second and third vectors.
  atomStruct atomA(1, 1, 0, 0);
  atomStruct atomB(1, 0, 1, 0);
  atomStruct atomC(1, 0, 0, 1);

  atomA = getAtomInCartCoords(atomA);
  atomB = getAtomInCartCoords(atomB);
  atomC = getAtomInCartCoords(atomC);

  vector<vector<double>> ret;

  double vecs[3][3];
  vecs[0][0] = atomA.x; vecs[0][1] = atomA.y; vecs[0][2] = atomA.z;
  vecs[1][0] = atomB.x; vecs[1][1] = atomB.y; vecs[1][2] = atomB.z;
  vecs[2][0] = atomC.x; vecs[2][1] = atomC.y; vecs[2][2] = atomC.z;

  // Sometimes it will express a small number as 1e-8 instead of 0
  // Express it as 0 instead
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      if (fabs(vecs[i][j]) < 1e-7) vecs[i][j] = 0;
    }
  }

  for (size_t i = 0; i < 3; i++)
    ret.push_back(vector<double>(vecs[i], vecs[i] + 3));

  return ret;
}

double Crystal::getVolume() const
{
  return m_lattice.a * m_lattice.b * m_lattice.c * getUnitVolume();
}

void Crystal::rescaleVolume(double newVolume)
{
  if (newVolume < 0) {
    cout << "Error! Crystal::rescaleVolume() was called to rescale the volume "
         << "to be a negative number, " << newVolume << "! Volume will not be "
         << "rescaled.\n";
    return;
  }

  double scalingFactor = pow(newVolume / getVolume(), 1.0/3.0);

  // Since atoms are all in fractional coordinates, we don't have to worry
  // about their positions changing...
  m_lattice.a *= scalingFactor;
  m_lattice.b *= scalingFactor;
  m_lattice.c *= scalingFactor;
  // This may put the lattice parameters outside the max and min bounds.
  // Might want to check up on those
}

double Crystal::getDistance(const atomStruct& as1,
                            const atomStruct& as2) const
{
  atomStruct cAs1 = getAtomInCartCoords(as1);
  atomStruct cAs2 = getAtomInCartCoords(as2);

  return sqrt(pow(cAs1.x - cAs2.x, 2.0) +
              pow(cAs1.y - cAs2.y, 2.0) +
              pow(cAs1.z - cAs2.z, 2.0));
}

double Crystal::findNearestNeighborAtomAndDistance(const atomStruct& as,
                                                   atomStruct& neighbor) const
{
  int ind = getAtomIndexNum(as);

  if (ind == -1) {
    cout << "Error: " << __FUNCTION__ << " was called for an atom that is not "
         << "a member of this cell!\n";
    return 0;
  }

  Crystal tempCrystal = *this;

  // We need to center the cell around this atom so that we don't run into the
  // problem of missing short distances caused by periodicity
  tempCrystal.centerCellAroundAtom(ind);

  size_t neighborInd = 0;
  double smallestDistance = 1000000.00;
  vector<atomStruct> tempAtoms = tempCrystal.getAtoms();
  for (size_t i = 0; i < tempAtoms.size(); i++) {
    if (i == ind) continue;
    double newDistance = getDistance(tempAtoms.at(ind), tempAtoms.at(i));
    if (newDistance < smallestDistance) {
      smallestDistance = newDistance;
      neighborInd = i;
    }
  }

  // Set the neighbor
  neighbor = m_atoms.at(neighborInd);

#ifdef NEAREST_NEIGHBOR_DEBUG
  cout << "Nearest neighbor is:\n";
  printAtomInfo(neighbor);

  cout << "distance is: " << smallestDistance << "\n";
#endif

  return smallestDistance;
}

int Crystal::getAtomIndexNum(const atomStruct& as) const
{
  for (size_t i = 0; i < m_atoms.size(); i++) {
    if (m_atoms.at(i) == as) return i;
  }
  cout << "Error in " << __FUNCTION__ << ": atom not found!\n";
  return -1;
}

void Crystal::centerCellAroundAtom(const atomStruct& as)
{
  int ind = getAtomIndexNum(as);

  if (ind == -1) {
    cout << "Error in " << __FUNCTION__ << ": a request was made to fill a "
         << "cell with an atom that is not a member of the cell!\n";
    return;
  }

  centerCellAroundAtom(ind);
}

void Crystal::centerCellAroundAtom(size_t ind)
{
  atomStruct& as = m_atoms.at(ind);

#ifdef CENTER_CELL_DEBUG
  cout << "Atom to be centered:\n";
  printAtomInfo(as);
  cout << "Before centering:\n";
  printAtomInfo();
#endif

  // Let's find the distances which we must shift the atoms so the one
  // at ind can be centered
  double dx, dy, dz;
  dx = 0.5 - as.x;
  dy = 0.5 - as.y;
  dz = 0.5 - as.z;

  for (size_t i = 0; i < m_atoms.size(); i++) {
    m_atoms[i].x += dx;
    m_atoms[i].y += dy;
    m_atoms[i].z += dz;
  }

  wrapAtomsToCell();

#ifdef CENTER_CELL_DEBUG
  cout << "After centering:\n";
  printAtomInfo();
#endif
}

// Radii should have already been scaled and set before calling this
double Crystal::getMinIAD(const atomStruct& as1, const atomStruct& as2) const
{
  double rad1 = ElemInfo::getRadius(as1.atomicNum, m_usingVdwRadii);
  double rad2 = ElemInfo::getRadius(as2.atomicNum, m_usingVdwRadii);
  return rad1 + rad2;
}

bool Crystal::areIADsOkay() const
{
  // We don't have to check the last atom if we checked all others
  for (size_t i = 0; i < m_atoms.size() - 1; i++) {
    if (!areIADsOkay(m_atoms.at(i))) return false;
  }
  return true;
}

bool Crystal::areIADsOkay(const atomStruct& as) const
{
  Crystal tempCrystal = *this;

  size_t ind = getAtomIndexNum(as);

  // We need to center the cell around this atom so that we don't run into the
  // problem of missing short distances caused by periodicity
  tempCrystal.centerCellAroundAtom(ind);

  vector<atomStruct> temp = tempCrystal.getAtoms();
  for (size_t i = 0; i < temp.size(); i++) {
    if (i == ind) continue;
    double minIAD = getMinIAD(temp.at(ind), temp.at(i));
    double dist = getDistance(temp.at(ind), temp.at(i));

    if (dist < minIAD) {
#ifdef IAD_DEBUG
      cout << "In " << __FUNCTION__ << ", minIAD failed!\n";
      cout << "  The distance is " << dist << " and the minIAD is " << minIAD
           << "\n";
      cout << "  Atoms responsible for failure are as follows:\n";
      printAtomInfo(as);
      printAtomInfo(neighbor);
#endif
      return false;
    }
  }
  return true;
}

// number of atoms and atomic number
typedef std::pair<uint, uint> numAndType;
vector<numAndType> getNumOfEachType(const vector<uint>& atoms)
{
  START_FT;
  vector<uint> atomsAlreadyCounted;
  vector<numAndType> numOfEachType;
  for (size_t i = 0; i < atoms.size(); i++) {
    size_t size = 0;
    // If we already counted this one, just continue
    if (std::find(atomsAlreadyCounted.begin(), atomsAlreadyCounted.end(),
                  atoms.at(i)) != atomsAlreadyCounted.end()) continue;
    for (size_t j = 0; j < atoms.size(); j++)
      if (atoms.at(j) == atoms.at(i)) size++;
    numOfEachType.push_back(make_pair(size, atoms.at(i)));
    atomsAlreadyCounted.push_back(atoms.at(i));
  }
  // Sort from largest to smallest
  sort(numOfEachType.begin(), numOfEachType.end(), greaterThan);
  return numOfEachType;
}

/* POSCAR format goes as such:
 *
 * Title
 * Scaling factor
 * Lattice vector for a
 * Lattice vector for b
 * Lattice vector for c
 * Element Symbols
 * Number of each element
 * Cartesian or Direct
 * Atom coordinates
 */
string Crystal::getPOSCARString(const string& title) const
{
  stringstream ss;

  // Set up the needed info
  ss << fixed << setprecision(15);
  vector<vector<double>> latticeVecs = getLatticeVecs();
  vector<numAndType> atomCounts = getNumOfEachType(getVectorOfAtomicNums());
  vector<string> symbols;
  for (size_t i = 0; i < atomCounts.size(); i++) {
    symbols.push_back(ElemInfo::getAtomicSymbol(atomCounts.at(i).second));
  }

  // Write to the POSCAR!
  ss << title << "\n"; // Title
  ss << "1.00000\n"; // Scaling factor

  for (size_t i = 0; i < 3; i++) {  // Lattice vectors
    for (size_t j = 0; j < 3; j++) {
      ss << " " << setw(20) << latticeVecs[i][j];
    }
    ss << "\n";
  }

  for (size_t i = 0; i < symbols.size(); i++) { // Symbols
    ss << "  " << setw(3) << symbols.at(i);
  }
  ss << "\n";

  for (size_t i = 0; i < atomCounts.size(); i++) { // Atom counts
    ss << "  " << setw(3) << atomCounts.at(i).first;
  }
  ss << "\n";

  ss << "Direct\n"; // We're just going to use fractional coordinates

  for (size_t i = 0; i < m_atoms.size(); i++) { // Atom coords
    ss << "  " << m_atoms.at(i).x << "  " << m_atoms.at(i).y << "  "
      << m_atoms.at(i).z << "\n";
  }

  // We're done!
  return ss.str();
}

void Crystal::writePOSCAR(const string& filename, const string& title) const
{
  ofstream f;
  f.open(filename);
  if (!f.is_open()) {
    cout << "Error in " << __FUNCTION__ << ": failed to open '"
         << filename << "' for writing!\n";
    return;
  }

  f << getPOSCARString(title);

  f.close();
}

bool Crystal::readPOSCAR(const string& filename,
                         latticeStruct& lattice,
                         vector<atomStruct>& atoms)
{
  ifstream f;
  f.open(filename);

  if (!f.is_open()) {
    cerr << "Error in " << __FUNCTION__ << ": file '" << filename
         << "' was not opened successfully! Please "
         << "check the permissions and that it exists.\n";
    return false;
  }

  std::string line;
  bool cart = false;

  // First line is comment
  getline(f, line);

  // Next line is scale factor
  getline(f, line);
  float scalingFactor = stof(line);

  float cell[3][3];

  // Next comes the matrix
  for (size_t i = 0; i < 3; i++) {
    getline(f, line);
    vector<string> theSplit = reduceAndSplit(line);
    // If this is not three, then there is some kind of error in the line
    if (theSplit.size() != 3) {
      cerr << "theSplit.size() is " << theSplit.size() << "\n";
      cerr << "In " << __FUNCTION__ << " while reading file '"
           << filename << "', error reading "
           << "the following line: '" << line << "'\n";
      return false;
    }
    cell[i][0] = (stof(theSplit.at(0)) * scalingFactor);
    cell[i][1] = (stof(theSplit.at(1)) * scalingFactor);
    cell[i][2] = (stof(theSplit.at(2)) * scalingFactor);
  }

  // Let's get a, b, c, alpha, beta, and gamma
  latticeStruct ls;
  ls.a = sqrt(pow(cell[0][0], 2) + pow(cell[0][1], 2) + pow(cell[0][2], 2));
  ls.b = sqrt(pow(cell[1][0], 2) + pow(cell[1][1], 2) + pow(cell[1][2], 2));
  ls.c = sqrt(pow(cell[2][0], 2) + pow(cell[2][1], 2) + pow(cell[2][2], 2));

  ls.alpha = rad2deg(acos(dotProduct<float>(cell[1], cell[2], 3) / (ls.b * ls.c)));
  ls.beta  = rad2deg(acos(dotProduct<float>(cell[0], cell[2], 3) / (ls.a * ls.c)));
  ls.gamma = rad2deg(acos(dotProduct<float>(cell[0], cell[1], 3) / (ls.a * ls.b)));

  // Sometimes, atomic symbols go here.
  getline(f, line);
  vector<uint> atomicNumbers;
  // If it contains atomic symbols, store them
  if (!containsOnlyIntsAndSpaces(line)) {
    vector<string> symbolsList = reduceAndSplit(line);
    for (size_t i = 0; i < symbolsList.size(); i++)
      atomicNumbers.push_back(ElemInfo::getAtomicNum(symbolsList.at(i)));
    // This next one should be atom types
    getline(f, line);
  }

  vector<string> theSplit = reduceAndSplit(line);
  vector<uint> atomCounts;
  for (size_t i = 0; i < theSplit.size(); i++) {
    uint atomicNum = stoi(theSplit.at(i));
    atomCounts.push_back(atomicNum);
  }

  // If we never filled up the atomic numbers, fill them up
  // now with "1, 2, 3..."
  if (atomicNumbers.size() == 0) {
    for (size_t i = 1; i <= atomCounts.size(); i++) atomicNumbers.push_back(i);
  }

  if (atomicNumbers.size() != atomCounts.size()) {
    cerr << "Error in " << __FUNCTION__ << ": the atom counts and atom "
         << "symbols do not match in number for POSCAR '" << filename
         << "'\n";
    return false;
  }

  // Starts with either [Ss]elective dynamics, [KkCc]artesian, or
  // other for fractional coords.
  getline(f, line);
  // If selective dynamics, get the next line
  if (line.at(0) == 'S' || line.at(0) == 's')
    getline(f, line);

  // Check if we're using cartesian or fractional coordinates:
  if (line.at(0) == 'K' || line.at(0) == 'k' ||
      line.at(0) == 'C' || line.at(0) == 'c' )
    cart = true;
  else
    cart = false;

  vector<atomStruct> tmpAtoms;
  for (size_t i = 0; i < atomCounts.size(); i++) {
    uint atomicNum = atomicNumbers.at(i);
    for (size_t j = 0; j < atomCounts.at(i); j++) {
      getline(f, line);
      vector<string> theSplit = reduceAndSplit(line);
      if (theSplit.size() != 3) {
        cerr << "Error in " << __FUNCTION__ << ": atomic coordinates are "
             << "not the right size in file '" << filename
             << "'\n";
        return false;
      }
      atomStruct tmp(atomicNum,
                     stof(theSplit.at(0)),
                     stof(theSplit.at(1)),
                     stof(theSplit.at(2)));
      // If we have cartesian coordinates, we need to adjust the positions
      // It's gonna look ugly...
      if (cart) {
        // We need alpha, beta, and gamma in radians
        double alpha = deg2rad(ls.alpha);
        double beta = deg2rad(ls.beta);
        double gamma = deg2rad(ls.gamma);
        // unit volume
        double v = pow(1.0 - pow(cos(alpha), 2.0) - pow(cos(beta), 2.0) - pow(cos(gamma), 2.0) +
                       2 * cos(alpha) * cos(beta) * cos(gamma), 0.5);
        tmp.x = (1.0 / ls.a) * tmp.x + (-cos(gamma) / (ls.a * sin(gamma))) * tmp.y +
                ((cos(alpha) * cos(gamma) - cos(beta)) / (ls.a * v * sin(gamma))) * tmp.z;
        tmp.y = (1.0 / (ls.b * sin(gamma))) * tmp.y +
                     ((cos(beta) * cos(gamma) - cos(alpha))
                      / (ls.b * v * sin(gamma))) * tmp.z;
        tmp.z = (sin(gamma) / (ls.c * v)) * tmp.z;
      }
      tmpAtoms.push_back(tmp);
    }
  }

  f.close();
  lattice = ls;
  atoms = tmpAtoms;
  return true;
}

string Crystal::getAtomInfoString(const atomStruct& as)
{
  stringstream s;
  s << "  For " << as.atomicNum << ", coords are: ("
    << as.x << "," << as.y << "," << as.z << ")\n";
  return s.str();
}

void Crystal::printAtomInfo(const atomStruct& as)
{
  cout << getAtomInfoString(as);
}

string Crystal::getAtomInfoString() const
{
  stringstream s;
  for (size_t i = 0; i < m_atoms.size(); i++) s << getAtomInfoString(m_atoms.at(i));
  return s.str();
}

void Crystal::printAtomInfo() const
{
  cout << getAtomInfoString();
}

string Crystal::getLatticeInfoString() const
{
  stringstream s;
  s << "a: " << m_lattice.a << "\n";
  s << "b: " << m_lattice.b << "\n";
  s << "c: " << m_lattice.c << "\n";
  s << "alpha: " << m_lattice.alpha << "\n";
  s << "beta: " << m_lattice.beta << "\n";
  s << "gamma: " << m_lattice.gamma << "\n";
  return s.str();
}

void Crystal::printLatticeInfo() const
{
  cout << getLatticeInfoString();
}

void Crystal::printLatticeVecs() const
{
  vector<vector<double>> vecs = getLatticeVecs();
  cout << "Lattice vecs are the following:\n";
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      cout << vecs[i][j] << " ";
    }
    cout << "\n";
  }
}

string Crystal::getCrystalInfoString() const
{
  stringstream s;
  s << "\n**** Printing Crystal Info ****\n";
  s << getLatticeInfoString() << getAtomInfoString();
  s << "**** End Crystal Info ****\n\n";
  return s.str();
}

void Crystal::printCrystalInfo() const
{
  cout << getCrystalInfoString();
}

void Crystal::printIADs() const
{
  vector<atomStruct> atoms = getAtoms();
  for (size_t i = 0; i < atoms.size(); i++) {
    cout << "For atom with index " << i << " and atomicNum " << atoms.at(i).atomicNum << ", the following are the neighbors:\n";
    Crystal tempCrystal = *this;

    // We need to center the cell around this atom so that we don't run into the
    // problem of missing short distances caused by periodicity
    tempCrystal.centerCellAroundAtom(i);

    size_t neighborInd = 0;
    double smallestDistance = 1000000.00;
    vector<atomStruct> tempAtoms = tempCrystal.getAtoms();
    for (size_t j = i + 1; j < tempAtoms.size(); j++) {
      double newDistance = getDistance(tempAtoms.at(i), tempAtoms.at(j));
      cout << "index " << j << " and atomicNum " << tempAtoms.at(j).atomicNum << ": " << newDistance << "\n";
    }
  }
}
