const int numSite = 10;
const int numEval = 10;
const double zMag = 0.0;
const std::string flagBoun = "PBC";

std::vector<int> bas;
const int dim = BasisConstruct(bas);

/*
BasisConstruct() function is to construct spin basis with conserved total Sz and return the dimension of the Hilbert space.
*/
int BasisConstruct(std::vector<int>& vec) {
        int max = 2 << (numSite-1);
        int dim = 0;
        for (int i = 0; i < max; i++) {
            std::bitset<numSite> b(i);
            if (b.count() == (unsigned int)(numSite / 2.0)) { // Impose the symmetry of total Sz = 0. 
                vec.push_back(i);
                dim++;
                }
            }
        return dim;
        }
