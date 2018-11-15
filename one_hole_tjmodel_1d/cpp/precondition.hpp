const int numSite = 10;
const int numHole = 1;
const int numEval = 10; // Desired eigenvalues to be found by ARPACK++.
const int numSam = 1;
const double step = 0.1;
const double flux = 0.0*PI;
const double zMag = 0.5;
const std::string flagBoun = "PBC";
const bool sigma = false; // Used to switch between t-J model and \sigma-t-J model.

std::vector<int> spinBasis;
std::vector<int> holeBasis;
const int subDim = SpinBasisConstruct(spinBasis);
const int holeDim = HoleBasisConstruct(holeBasis); 
const int dim = subDim*holeDim;

/*
BasisConstruct() function is to construct spin basis with conserved total Sz and return the dimension of the Hilbert space.
*/
int SpinBasisConstruct(std::vector<int> &vec) {
        int max = int(pow(2, numSite-1)); // Possible maxmium dimension of the Hilbert space for the sub-space in terms of pure spin configuration with the hole fixed.
        int dim = 0;
        for (int i = 0; i < max; ++i) {
            std::bitset<numSite> b(i);
            if (b.count() == (unsigned int)(zMag+(numSite-numHole)/2.0)) {
                vec.push_back(i);
                dim++;
                }
            }
        return dim;
        }

int HoleBasisConstruct(std::vector<int> &vec) {
        int dim = 0;
        for (int i = 0; i < numSite; ++i) {
            vec.push_back(i);
            dim++;
            }
        return dim;
        }
