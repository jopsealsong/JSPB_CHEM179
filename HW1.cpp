//
//  main.cpp
//  CHEM179_HW1
//
//  Created by Joseph Song on 2/2/25.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <sstream>
#include <vector>
using namespace std;

/* ----- PART 1: PARSING INPUT FILE ----- */

/* stores input parameters in class */
struct Atom {
public:
    double x,y,z;
    
    Atom(double x, double y, double z) : x(x), y(y), z(z) {}
};

/* parses file to create vector of Atoms. Used AI for syntax guidance and
 troubleshooting errors */
std::vector<Atom> initiateAtoms(const std::string& filePath) {
    std::vector<Atom> atoms;
    std::ifstream file(filePath);
    
    if (!file) {
        std::cerr << "Cannot open file" << std::endl;
        return atoms;
    }
    
    std::string line;
    int numberAtoms;
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        if (!(ss >> numberAtoms) || numberAtoms <= 0) {
            std::cerr << "ERROR: INVALID ATOM COUNT" << std::endl;
            return atoms;
        }
        atoms.resize(numberAtoms);
    } else {
        std::cerr << "ERROR: EMPTY FILE" << std::endl;
        return atoms;
    }
    
    int lineNumber = 2;
    int actualLines = 0;
    
    while (std::getline(file, line)) {
        actualLines++;
        
        std::stringstream ss(line);
        int atomicNumber;
        double x, y, z;
        
        if (!(ss >> atomicNumber >> x >> y >> z)) {
            std::cerr << "ERROR: IMPROPER FORMATTING" << std::endl;
            return {};
        }
        if (atomicNumber != 79) {
            std::cerr << "ERROR: NOT GOLD ATOMS" << std::endl;
            return {};
        }
        
        Atom atom = Atom(x, y, z);
        atoms.push_back(atom);
        lineNumber++;
    }
    
    file.close();
    
    if (actualLines != numberAtoms) {
        std::cerr << "ERROR: IMPROPER NUMBER OF ATOMS" << std::endl;
    }
    
    return atoms;
}
    

/* ----- PART 2: ENERGY CALCULATIONS ----- */

/* computes distance between two atoms */
double computeDistance(Atom atom1, Atom atom2)
{
    return sqrt(pow(atom1.x - atom2.x, 2) + pow(atom1.y - atom2.y, 2) + pow(atom1.z - atom2.z, 2));
}

/* computes Lennard Jones energy between two atoms */
double lennardJones(double R)
{
    double sigma = 2.951;
    double epsilon = 5.290;
    double LJterm = sigma/R;
    double energy = epsilon*(pow(LJterm, 12) - 2*pow(LJterm, 6));
    return energy;
};

/* iterates through pairs of atoms to compute total energy of a given cluster config */
double computeEnergy(std::list<Atom> atoms)
{
    double energy = 0.000;
    int curr_atom = 0;
    for (std::list<Atom>::iterator atomIterator = atoms.begin(); atomIterator != atoms.end(); atomIterator++)
    {
        std::list<Atom>::iterator atomPair = atoms.begin();
        for (int i = 0; i <= curr_atom; i++) {
            atomPair++;
        }
        while (atomPair != atoms.end()) {
            energy += lennardJones(computeDistance(*atomIterator, *atomPair));
            atomPair++;
        }
        curr_atom++;
    }
    return energy;
}

/* ------ PART 3: FORCE CALCULATIONS AND FINITE DIFFERENCES ----- */

/* computes analytical derivative of LJ */
double analyticalDeriv(Atom atom1, Atom atom2, double Atom::*variable)
{
    double sigma = 2.951;
    double epsilon = 5.290;
    double R = computeDistance(atom1, atom2);
    double firstTerm = epsilon*(12*(pow(sigma, 12)/pow(R, 13)) - 12*(pow(sigma, 6)/pow(R, 7)));
    double secondTerm = (atom2.*variable - atom1.*variable)/R;
    return firstTerm*secondTerm;
}

/* computes force applied on an atom in a cardinal direction */
double computeAnalyForce(Atom atom, std::list<Atom> atoms, double Atom::*dim)
{
    double force = 0.000;
    for (Atom& curr_atom : atoms) {
        if (&curr_atom != &atom) {
            force += analyticalDeriv(atom, curr_atom, dim);
        }
    }
    return force;
}

/* approximates force via finite forward difference */
double computeFiniteForward(double step, std::list<Atom> atoms, int index, double Atom::*variable)
{
    double initialEnergy = computeEnergy(atoms);
    std::list<Atom> atoms1 = atoms;
    std::list<Atom>::iterator atomPair = atoms1.begin();
    for (int i = 0; i < index; i++)
    {
        atomPair++;
    }
    *atomPair.*variable += step;
    double finalEnergy = computeEnergy(atoms1);
    return -(finalEnergy - initialEnergy)/step;
}

/* approximates force via finite central difference */
double computeCentralForward(double step, std::list<Atom> atoms, int index, double Atom::*variable)
{
    std::list<Atom> atoms1 = atoms;
    std::list<Atom>::iterator atomPair = atoms1.begin();
    for (int i = 0; i < index; i++)
    {
        atomPair++;
    }
    *atomPair.*variable += step;
    double forwardEnergy = computeEnergy(atoms1);
    *atomPair.*variable -= 2*step;
    double backwardEnergy = computeEnergy(atoms1);
    return -(forwardEnergy - backwardEnergy)/(2*step);
}

/* ----- PART 4: STEEPEST DESCENT ----- */

/* computes gradient vector (aka force vector) */
std::vector<double> computeGradient(std::list<Atom>& atoms) {
    std::vector<double> dimensions;
    dimensions.reserve(atoms.size()*3);
    for (Atom& atom : atoms) {
        dimensions.push_back(computeAnalyForce(atom, atoms, &Atom::x));
        dimensions.push_back(computeAnalyForce(atom, atoms, &Atom::y));
        dimensions.push_back(computeAnalyForce(atom, atoms, &Atom::z));   
    }
    return dimensions;
};

/* computes gradient modulus */
double gradientLength(std::vector<double> gradient) {
    double sum = 0.000;
    for (double dim : gradient) {
        sum += pow(dim, 2);
    }
    return sqrt(sum);
};

/* adjusts gradient modulus by a inputted value */
std::vector<double> adjustGradMagnitude(std::vector<double> gradient, double value) {
    for (double& dim : gradient) {
        dim *= value;
    }
    return gradient;
}

/* takes a step toward a specified direction and magnitude in the PES */
void incrementPos(std::list<Atom>& atoms, std::vector<double> descent) {
    if (descent.size() != 3 * atoms.size()) {
        throw std::invalid_argument("array sizes don't match");
    }
    
    std::vector<double>::iterator dimIter = descent.begin();
    for (Atom& atom : atoms) {
        atom.x -= *dimIter;
        dimIter++;
        atom.y -= *dimIter;
        dimIter++;
        atom.z -= *dimIter;
        dimIter++;
    }
};

/* implementation of steepest descent algorithm without line search */
void steepestDescent(std::list<Atom>& atoms, const double& tol, double& step)
{
    std::vector<double> gradient = computeGradient(atoms);
    while (gradientLength(gradient) >= tol) {
        std::vector<double> stepdir = gradient;
        stepdir = adjustGradMagnitude(stepdir, step/gradientLength(stepdir));
        incrementPos(atoms, stepdir);
        
        std::vector<double> fgradient = computeGradient(atoms);
        if (gradientLength(fgradient) >= gradientLength(gradient)) {
            step *= 1.10;
        } else {
            step /= 2;
        }
        gradient = fgradient;
    }
};

/* ----- PART 5: STEEPEST DESCENT WITH LINE SEARCH ----- */

void goldsection(std::vector<double>& interval, std::list<Atom>& point_one, std::list<Atom>& point_two, double tol) {
    const double goldratio = 0.618;
    double intervalLength = gradientLength(interval);
    
    if (intervalLength <= tol) {
        return;
    }
    
    interval = adjustGradMagnitude(interval, goldratio);
        
    std::list<Atom> point_three = point_two;
    incrementPos(point_three, adjustGradMagnitude(interval, -1));
    std::list<Atom> point_four = point_one;
    incrementPos(point_four, interval);
        
    if (computeEnergy(point_four) > computeEnergy(point_three)) {
        point_two = point_four;
    } else {
        point_one = point_three;
    }
    goldsection(interval, point_one, point_two, tol);
}

void DescentLineSearch(std::list<Atom>& atoms, const double& tol, double& step)
{
    std::vector<double> gradient = computeGradient(atoms);
    while (gradientLength(gradient) >= tol) {
        std::vector<double> stepdir = gradient;
        stepdir = adjustGradMagnitude(stepdir, step/gradientLength(stepdir));
        
        std::list<Atom> point_one = atoms;
        double E1 = computeEnergy(point_one);
        
        std::list<Atom> point_two = point_one;
        incrementPos(point_two, adjustGradMagnitude(stepdir, 1/2));
        double E2 = computeEnergy(point_two);
        
        std::list<Atom> point_three = point_two;
        incrementPos(point_three, adjustGradMagnitude(stepdir, 1/2));
        double E3 = computeEnergy(point_three);
        
        if (E1 > E2 & E3 > E2) {
            goldsection(stepdir, point_one, point_two, tol);
            incrementPos(atoms, adjustGradMagnitude(stepdir, 1/2));
        } else {
            incrementPos(atoms, stepdir);
        }
        
        std::vector<double> fgradient = computeGradient(atoms);
        if (gradientLength(fgradient) >= gradientLength(gradient)) {
            step *= 1.10;
        } else {
            step /= 2;
        }
        gradient = fgradient;
    }
};

int main(int argc, const char * argv[]) {
    
    std::cout << "Hello, World!\n";
    return 0;
}
