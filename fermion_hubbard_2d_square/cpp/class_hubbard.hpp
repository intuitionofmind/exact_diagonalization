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
        void SaveHilbert();
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
                        ++dim;
                        }
                    }
                }
            }
        std::sort(mBasis.begin(), mBasis.end());
        mDim = dim;
        
/*        for (int i = 0; i < dim; ++i) {
            int s = mBasis[i];
            std::cout << s << " ";
            boost::dynamic_bitset<> bit(numSite*2, s);
            for (int l = 0; l < numSite*2; ++l) { std::cout << bit[l]; }
            std::cout << std::endl;
            }*/
        }

// Class destructor.
template<typename T>
Hubbard<T>::~Hubbard() {}

template<typename T>
int Hubbard<T>::HilbertDim() { return mDim; }

template<typename T>
void Hubbard<T>::SaveHilbert() {
        std::ofstream file_hilbert("hilbert_space", std::ios_base::app);
        for (int l = 0; l < mDim; ++l) { file_hilbert << mBasis[l] << ","; }
        file_hilbert.close();
        }

template<typename T>
int Hubbard<T>::Coordinate(int x, int y) { return y*mNumSiteX+x; }

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
                            T fSign = pow(-1.0, cross);
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
                                T fSign = pow(-1.0, cross);
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
                                T fSign = pow(-1.0, cross);
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

// Note that the ArcomStdEig() in ARPACKPP will not sort the eigenvalues you want while ArsymStdEig() does. SortEval() funtion helps you to sort the eigenvalues and its order is stored in the vector "order" such that you can access the i'th smallest eigenvalues according to "order[i]" in the original sequence. 
/*
template<typename T>
void Hubbard<T>::Sort(int n, T* w, T* v) { 
        std::vector<double> vec;
        for (int i = 0; i < n; ++i) { vec.push_back(std::real(w[i])); }
        std::sort(vec.begin(), vec.end());

        auto temp = new T[mDim*n];

        for (int i = 0; i < n; ++i) {
            std::vector<int>::iterator it = std::find(vec.begin(), vec.end(), std::real(w[i])); 
            w[i] = T(vec[i]);
            }


        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::vector<int>::iterator it = std::find(order.begin(), order.end(), j); 
                if (std::real(w2[i]) == std::real(w1[j]) && it == order.end()) {
                    order.push_back(j);
                    }
                }
            }

        delete [] temp;
        }
    */

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
                std::cout << std::real(Dot(v, u)) << ", ";
                }
            std::cout << std::endl;
            }
        std::cout << std::endl;
        delete [] u;
        delete [] v;
        delete [] w;
        }
