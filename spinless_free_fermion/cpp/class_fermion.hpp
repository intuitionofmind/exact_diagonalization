#include "head.hpp"

// A class is an object. Here the object FreeFermion denotes a quantum system consisting of a specific lattce and kind of Hamiltoian.
template<typename T>
class FreeFermion {
        private:
        // Paramaters.
        int mNum; // Number of electrons. 
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
        FreeFermion (int num, int numSiteX, int numSiteY, std::string bouConX, std::string bouConY); // Constructor.
        ~FreeFermion (); // Destructor.

        int HilbertDim();
        int Coordinate(int x, int y);
        void Hamiltonian(T* v, T* w);
        
        void SetOne(T* v, int i);
        T Dot(T* v, T* w);
        void PrintHam();
        };

template<typename T>
FreeFermion<T>::FreeFermion (int num, int numSiteX, int numSiteY, std::string bouConX, std::string bouConY) {
        mNum = num;
        mNumSiteX = numSiteX;
        mNumSiteY = numSiteY;
        mBouConX = bouConX;
        mBouConY = bouConY;

        // Construct the basis of the Hilbert space constrained by the U(1) symmetry.
        int numSite = numSiteX*numSiteY;
        long int max = int(pow(2, numSite));
        int dim = 0;
        for (int i = 0; i < max; ++i) {
            boost::dynamic_bitset<> bit(numSite, i);
            if (int(bit.count()) == num) {
                mBasis.push_back(bit.to_ulong());
                ++dim;
                }
            }
        std::sort(mBasis.begin(), mBasis.end());
        mDim = dim;
        std::cout << mDim << std::endl;
        }

// Class destructor.
template<typename T>
FreeFermion<T>::~FreeFermion() {}

template<typename T>
int FreeFermion<T>::Coordinate(int x, int y) { return y*mNumSiteX+x; }

template<typename T>
int FreeFermion<T>::HilbertDim() { return mDim; }

// Required by Arpack++ package handbook: There only requirements make by ARPACK++ are that member funtion Hamiltonian() musth have two pointers to vectors of type T as paraments and the input vector must precede the output vector.
template<typename T>
void FreeFermion<T>::Hamiltonian(T* v, T* w) {
        int numSite = mNumSiteX*mNumSiteY;

        for (int l = 0; l < mDim; ++l) { w[l] = 0.0; }
        for (int l = 0; l < mDim; ++l) {
            if (0.0 == v[l]) { continue; }
            int b = mBasis[l];
            boost::dynamic_bitset<> config(numSite, b);

            // Hamiltonian operation loop through all sites.
            for (int j = 0; j < mNumSiteY; ++j) {
                for (int i = 0; i < mNumSiteX; ++i) {

                    for (int c = 0; c < 2; ++c) {
                        bool flagX = false; // Mark the boundary terms.
                        bool flagY = false;
                        int current = Coordinate(i, j); // Location of the spin in the representation string.
                        int forward = 0;
                        switch (c) {
                            case 0:
                            forward = Coordinate((i+1) % mNumSiteX, j);
                            if (forward < current) { flagX = true; }
                            break;
                            case 1:
                            forward = Coordinate(i, (j+1) % mNumSiteY);
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
                            std::vector<int>::iterator it = std::lower_bound(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));

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
                                std::vector<int>::iterator it = std::lower_bound(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));
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
                                std::vector<int>::iterator it = std::lower_bound(mBasis.begin(), mBasis.end(), int(temp.to_ulong()));
                                w[std::distance(mBasis.begin(), it)] += -1.0*fSign*v[l];
                                }
                            }
                        }
                    }
                }
            }
        }

template<typename T>
void FreeFermion<T>::SetOne(T* v, int i) {
        for (int l = 0; l < mDim; ++l) { v[l] = 0.0; }
        v[i] = 1.0;
        }

template<typename T>
T FreeFermion<T>::Dot(T* v, T* w) {
        T r = 0.0;
        for (int l = 0; l < mDim; ++l) { r += std::conj(v[l])*w[l]; }
        return r;
        }

template<typename T>
void FreeFermion<T>::PrintHam() {
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
