#include "head.hpp"

// A class is an object. Here the object Hubbard denotes a quantum system consisting of a specific lattce and kind of Hamiltoian.
template<typename T>
class Hubbard {
        private:
        // Paramaters.
        double mU; // Hubbard U. t=1.0 is set as energy unit. 
        int mNumUp; // Number of spin-up electrons. 
        int mNumDown; // Number of spin-down electrons.
        // Define the lattice.
        int mNumSiteX;
        int mNumSiteY;
        std::string mBouConX;
        std::string mBouConY;
        // Dimension and basis of the Hilbert.
        int mDim;
        std::vector<int> mBasis;

        public:
        Hubbard (double u, int numUp, int numDown, int numSiteX, int numSiteY, std::string bouConX, std::string bouConY); // Constructor.
        ~Hubbard (); // Destructor.

        int HilbertDim();
        int Coordinate(int x, int y);
        void Hamiltonian(T* v, T* w);
        
        void SetOne(T* v, int i);
        T Dot(T* v, T* w);
        void PrintHam();
        };

// To define a specific quantum system and build up the corresponding Hilbert space in the tensor product representation.
template<typename T>
Hubbard<T>::Hubbard (double u, int numUp, int numDown, int numSiteX, int numSiteY, std::string bouConX, std::string bouConY) {
        mU = u;
        mNumUp = numUp;
        mNumDown = numDown;
        mNumSiteX = numSiteX;
        mNumSiteY = numSiteY;
        mBouConX = bouConX;
        mBouConY = bouConY;

        // Construct the basis of the Hilbert space constrained by the U(1) symmetry.
        int numSite = numSiteX*numSiteY;
        long int max = int(pow(2, numSite));
        int dim = 0;
        for (int i = 0; i < max; ++i) {
            boost::dynamic_bitset<> bitUp(numSite, i);
            if (int(bitUp.count()) == numUp) {
                for (int j = 0; j < max; ++j) {
                    boost::dynamic_bitset<> bitDown(numSite, j);
                    if (int(bitDown.count()) == numDown) {
                        boost::dynamic_bitset<> bit;
                        for (int k = 0; k < numSite; ++k) { bit.push_back(bitUp[k]); } // Basis of tensor product of up- and down-spins.
                        for (int k = 0; k < numSite; ++k) { bit.push_back(bitDown[k]); }
                        mBasis.push_back(bit.to_ulong());
                        std::cout << bit.to_ulong() << " ";
                        for (int l = 0; l < numSite*2; ++l) { std::cout << bit[l]; }
                        std::cout << std::endl;
                        ++dim;
                        }
                    }
                }
            }
        mDim = dim;
        std::cout << mDim << std::endl;
        }

// Class destructor.
template<typename T>
Hubbard<T>::~Hubbard() {}

template<typename T>
int Hubbard<T>::Coordinate(int x, int y) { return y*mNumSiteX+x; }

template<typename T>
int Hubbard<T>::HilbertDim() { return mDim; }

// Required by Arpack++ package handbook: There only requirements make by ARPACK++ are that member funtion Hamiltonian() musth have two pointers to vectors of type T as paraments and the input vector must precede the output vector.
template<typename T>
void Hubbard<T>::Hamiltonian(T* v, T* w) {
        int numSite = mNumSiteX*mNumSiteY;

        for (int l = 0; l < mDim; ++l) { w[l] = 0.0; }
        for (int l = 0; l < mDim; ++l) {
            if (0.0 == v[l]) { continue; }
            int b = mBasis[l];
            boost::dynamic_bitset<> config(numSite*2, b);

            // Hamiltonian operation loop through all sites.
            for (int j = 0; j < mNumSiteY; ++j) {
                for (int i = 0; i < mNumSiteX; ++i) {
                    int k = Coordinate(i, j);
                    if (config[k] & config[numSite+k]) { w[l] += mU*v[l]; } // Diagonal Hubbard U onsite interaction.

                    for (int c = 0; c < 4; ++c) { // Four cases: 0: spin-up along x-direction; 1: spin-up along y-direction; 2: spin-down along x-direction; 3: spin-down along y-direction.
                        bool flagX = false; // Mark the boundary terms.
                        bool flagY = false;
                        int current = 0; // Location of the spin in the representation string.
                        int forward = 0;
                        switch (c) {
                            case 0:
                            current = k;
                            forward = Coordinate((i+1) % mNumSiteX, j);
                            if (forward < current) { flagX = true; }
                            break;
                            case 1:
                            current = k;
                            forward = Coordinate(i, (j+1) % mNumSiteY);
                            if (forward < current) { flagY = true; }
                            break;
                            case 2:
                            current = numSite+k;
                            forward = numSite+Coordinate((i+1) % mNumSiteX, j);
                            if (forward < current) { flagX = true; }
                            break;
                            case 3:
                            current = numSite+k;
                            forward = numSite+Coordinate(i, (j+1) % mNumSiteY);
                            if (forward < current) { flagY = true; }
                            break;
                            }

                        // For the hopping process within the graph, exclusive boundary hoppings. Note that it does not mean electrons hop from "current" to "forward" but the bond "current-forward".
                        if ((forward > current) && (config[current] ^ config[forward])) {
                            int cross = 0; // Count how many fermions crossed by the operator.
                            for (int m = current+1; m < forward; ++m) { if (config[m]) { ++cross; } }
                            double fSign = pow(-1.0, cross);
                            boost::dynamic_bitset<> temp(config);
                            temp.flip(current);
                            temp.flip(forward);
                            std::vector<int>::iterator it = std::find(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));
/*
                            for (int p = 0; p < numSite*2; ++p) { std::cout << config[p]; }
                            std::cout << std::endl;
                            for (int p = 0; p < numSite*2; ++p) { std::cout << temp[p]; }
                            std::cout << std::endl;
                            std::cout << config.to_ulong() << " " << temp.to_ulong() << " " << std::distance(mBasis.begin(), it) << " " << cross << std::endl;
                            */

                            w[std::distance(mBasis.begin(), it)] += -1.0*fSign*v[l];
                            }
                        // For the boundary hoppings.
                        else if (flagX && "PBC" == mBouConX) {
                            if (config[current] ^ config[forward]) {
                                int cross = 0;
                                for (int m = forward+1; m < current; ++m) { if (config[m]) { ++cross; } }
                                double fSign = pow(-1.0, cross);
                                boost::dynamic_bitset<> temp(config);
                                temp.flip(current);
                                temp.flip(forward);
                                std::vector<int>::iterator it = std::find(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));
                                w[std::distance(mBasis.begin(), it)] += -1.0*fSign*v[l];
                                }
                            }
                        else if (flagY && "PBC" == mBouConY) {
                            if (config[current] ^ config[forward]) {
                                int cross = 0;
                                for (int m = forward+1; m < current; ++m) { if (config[m]) { ++cross; } }
                                double fSign = pow(-1.0, cross);
                                boost::dynamic_bitset<> temp(config);
                                temp.flip(current);
                                temp.flip(forward);
                                std::vector<int>::iterator it = std::find(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));
                                w[std::distance(mBasis.begin(), it)] += -1.0*fSign*v[l];
                                }
                            }
                        }
                    }
                }
            }
        }

template<typename T>
void Hubbard<T>::SetOne(T* v, int i) {
        for (int l = 0; l < mDim; ++l) { v[l] = 0.0; }
        v[i] = 1.0;
        }

template<typename T>
T Hubbard<T>::Dot(T* v, T* w) {
        T r = 0.0;
        for (int l = 0; l < mDim; ++l) { r += std::conj(v[l])*w[l]; }
        return r;
        }

template<typename T>
void Hubbard<T>::PrintHam() {
        auto u = new T[mDim];
        auto v = new T[mDim];
        auto w = new T[mDim];
        for (int i = 0; i < mDim; ++i) {
            for (int j = 0; j < mDim; ++j) {
                SetOne(v, i);
                SetOne(w, j);
                Hamiltonian(w, u);
                std::cout << std::real(Dot(v, u)) << " "; 
                }
            std::cout << std::endl;
            }
        std::cout << std::endl;
        delete [] u;
        delete [] v;
        delete [] w;
        }
