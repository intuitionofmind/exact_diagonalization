/*
 * This it the file to define the Hamiltonian operoator and related operations.
 */

#include "arcomp.h"

template<int N>
int RetrieveSpin(std::bitset<N> config, int h) {
        std::bitset<N> spin;
        for (int i = 0; i < N; ++i) {
            if (i < h) { spin[i] = config[i]; }
            else { spin[i] = config[i+1]; }
            }
        return spin.to_ulong();
        }

// note that std::swap() function does not apply to std::bitset
template<int N>
void SwapBit(std::bitset<N>& config, int i, int j) {
        int a = config[i];
        config[i] = config[j];
        config[j] = a;
        }

template<typename T>
class tjChainHalf {
        private:
        int n;
        T J;
        std::vector<T> p;

        public:
        int Dim();
        T CoupStren();
        T Dot(T* v1, T* v2);
        void VecPlus(T* v1, T* v2, T* w);
        void VecMinus(T* v1, T* v2, T* w);
        void MultNumber(T a, T*v);
        void Copy(T* v1, T* v2);
        void Normalize(T* v);
        void SetZero();
        void SetOne(T* v, int i);
        void SortEval(int n, T* w1, T* w2, std::vector<int>& order);
        void MultVec(T* v, T* w);
        void TimeEvolution(double step, int num, T* vInit, T* vFinal);
        void Translation(T* v1, T* v2);
        T Correlation(T* v, int i, int j);
        void Marshall(T* v1, T* v2);

        // constructor
        tjChainHalf(int d, T j) { 
            n = d;
            J = j;
            std::vector<T> p(d);
            }  
        };

template<typename T>
int tjChainHalf<T>::Dim() { return n; }

template<typename T>
T tjChainHalf<T>::CoupStren() { return J; }
template<typename T> T tjChainHalf<T>::Dot(T* v1, T* v2) {
        // int len = Dim();
        T res = 0.0;
        // arcomplex<double> res (0.0, 0.0);
        // for (int i = 0; i < len; ++i) { res = res + std::conj(v1[i])*v2[i]; }
        return res;
        }

template<typename T>
void tjChainHalf<T>::VecPlus(T* v1, T* v2, T* w) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { w[i]= v1[i]+v2[i]; }
        }

template<typename T>
void tjChainHalf<T>::VecMinus(T* v1, T* v2, T* w) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { w[i]= v1[i]-v2[i]; }
        }

template<typename T>
void tjChainHalf<T>::MultNumber(T a, T* v) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v[i] = a*v[i]; }
        }

template<typename T>
void tjChainHalf<T>::Copy(T* v1, T* v2) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v2[i] = v1[i]; }
        }

template<typename T>
void tjChainHalf<T>::Normalize(T* v) {
        T z = Dot(v, v);
        double f = sqrt(std::real(z));
        int len = Dim();
        for (int i = 0; i < len; ++i) { v[i] = v[i]/f; }
        }

template<typename T>
void tjChainHalf<T>::SetZero() {
        for (auto val:p) { val = 0.; }
        // for (auto iter = p.begin(); iter != p.end(); ++iter) {}
        }

template<typename T>
void tjChainHalf<T>::SetOne(T* v, int i) {
        int len = Dim();
        for (int j = 0; j < len; ++j) { v[j] = 0.0; }
        v[i] = 1.0;
        }

/*
 * Note that the ArcomStdEig() in ARPACKPP will not sort the eigenvalues you want while ArsymStdEig() does. SortEval() funtion helps to sort the eigenvalues and its order is stored in the vector "order" such that you can access the i'th smallest eigenvalues according to "order[i]" in the original sequence. 
 */

template<typename T>
void tjChainHalf<T>::SortEval(int n, T* w1, T* w2, std::vector<int>& order) { 
        std::vector<double> vec;
        for (int i = 0; i < n; ++i) { vec.push_back(std::real(w1[i])); }
        std::sort(vec.begin(), vec.end());
        for (int i = 0; i < n; ++i) { w2[i] = T(vec[i]); }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::vector<int>::iterator it = std::find(order.begin(), order.end(), j); 
                if (std::real(w2[i]) == std::real(w1[j]) && it == order.end()) {
                    order.push_back(j);
                    }
                }
            }
        }

/*
 * One dimensional spin-1/2 one hole doped t-J Hamiltonian. Note that there is a 1/4 shift in the Hamiltonian: H = \sum_{\langle i, j\rangle}(S_{i}S_{j}-1/4)+t*c_{i}^{\dagger}c_{j}.
 * 
 * Required by Arpack++ package handbook:
 * There only requirements make by ARPACK++ are that MultVec musth have two pointers to vectors of type T as paraments and the input vector must precede the output vector.
 */

template<typename T>
void tjChainHalf<T>::MultVec(T* v, T* w) {
        int len = Dim();
        T J = CoupStren();
        for (int i = 0; i < len; ++i) { w[i] = 0.0; } 
        for (int i = 0; i < len; ++i) {
            if (0.0 == v[i]) { continue; }
            int s = spinBasis[i % subDim]; // convention is i = h*subDim+s  
            int h = holeBasis[int(i/subDim)];
            std::bitset<numSite> spinConfig(s); // the highest bits for the number of (numSite-numHole) are filled with '0' while do not have real meanings but just for convenience
            std::bitset<numSite> config; // initialize config with all 0s

            // Fuse hole and spin to give rise to a tJ configuration.
            int jj = 0;
            for (int j = 0; j < numSite; ++j) {
                if (h == j) { continue; }
                else {
                    config[j] = spinConfig[jj];
                    ++jj;
                    }
                }

            for (int j = 0; j < numSite; ++j) {
                double sign = 1.0; // sign for sigma tJ model
                T phase = 1.0; // sign for fermion hopping and possible flux on the boundary
                int jj = (j+1) % numSite;
                if (jj < j && "OBC" == flagBoun) { continue; }

                // Heisenberg J term
                if (j != h && jj != h && config[j] != config[jj]) {
                    w[i] -= 0.5*J*v[i]; // digonal Sz term 
                    std::bitset<numSite> temp (config);
                    temp.flip(j);
                    temp.flip(jj);
                    int ss = RetrieveSpin<numSite>(temp, h);
                    std::vector<int>::iterator k = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                    w[h*subDim+std::distance(spinBasis.begin(), k)] += 0.5*J*v[i];
                    }

                // hopping t term
                else if (jj == h) {
                    int hh = j;
                    if (sigma && 0 == config[j]) { sign = -1.0; }
                    // if (0 == jj) { phase = pow(-1.0, numSite)*std::polar(1.0, flux); } // across the boundary 
                    std::bitset<numSite> temp (config);
                    SwapBit<numSite>(temp, j, jj);
                    int ss = RetrieveSpin<numSite>(temp, hh);
                    std::vector<int>::iterator k = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                    w[hh*subDim+std::distance(spinBasis.begin(), k)] -= 1.0*sign*phase*v[i];
                    }
                else if (j == h) {
                    int hh = jj;
                    if (sigma && 0 == config[jj]) { sign = -1.0; }
                    // if (0 == jj) { phase = pow(-1.0, numSite)*std::polar(1.0, flux); }
                    std::bitset<numSite> temp (config);
                    SwapBit<numSite>(temp, j, jj);
                    int ss = RetrieveSpin<numSite>(temp, hh);
                    std::vector<int>::iterator k = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
                    w[hh*subDim+std::distance(spinBasis.begin(), k)] -= 1.0*sign*phase*v[i];
                    }   
                }
            }
        }

// Integrate the time evolution operator U=exp(-i*H*t) by num steps.
template<typename T>
void tjChainHalf<T>::TimeEvolution(double step, int num, T* vInit, T* vFinal) {
        arcomplex<double> imagUnit (0.0, 1.0);
        // In askar-Goldberg method, u, v, w denote n-1, n, n+1 states.
        auto u = new arcomplex<double>[dim];
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        auto temp = new arcomplex<double>[dim];
        auto tempTemp = new arcomplex<double>[dim];

        for (int i = 0; i < num; ++i) {
            if (0 == i) {
                Copy(vInit, w);
                Copy(w, v);
                }
            else if (1 == i) {
                MultVec(v, temp);
                MultNumber(imagUnit*step, temp);
                VecMinus(v, temp, w);
                Normalize(w);
                // Move a step forward.
                Copy(v, u);
                Copy(w, v);
                }
            else {
                MultVec(v, temp);
                MultVec(temp, tempTemp);
                MultNumber(2.0*imagUnit*step, temp);
                VecMinus(u, temp, w);
                MultVec(tempTemp, temp);  // For H^{3}. 
                MultNumber(0.5*imagUnit*step*step*step, temp);
                VecPlus(w, temp, w);
                Copy(v, u);
                Copy(w, v);
                }
            }

        Copy(w, vFinal);
        delete [] u;
        delete [] v;
        delete [] w;
        delete [] temp;
        delete [] tempTemp;
        }

template<typename T>
void tjChainHalf<T>::Translation(T* v1, T* v2) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { v2[i] = 0.0; }
        for (int i = 0; i < len; ++i) {
            int s = spinBasis[i % subDim]; // The convention is i = h*subDim+s.
            int h = holeBasis[int(i/subDim)];
            std::bitset<numSite> spinConfig(s);
            std::bitset<numSite> configTemp(s);
            configTemp.reset();
            int hh = (h+1) % numSite;
            double sign;
            int ss;
            if (h == numSite-1) { // This is a specific case in which the spin configuration is not changed
                sign = 1.0;
                ss = s;
                }
            else {
                sign = pow(-1.0, numSite);
                for (int j = 1; j < numSite-1; ++j) { configTemp[j] = spinConfig[j-1]; }
                configTemp[0] = spinConfig[numSite-2];
                ss = configTemp.to_ulong();
                }
            std::vector<int>::iterator k = std::lower_bound(spinBasis.begin(), spinBasis.end(), ss);
            v2[hh*subDim+std::distance(spinBasis.begin(), k)] = sign*v1[i];
            }
        }

template<typename T>
T tjChainHalf<T>::Correlation(T* v, int i, int j) {
        int len = Dim();
        T res = 0.0;
        for (int k = 0; k < len; ++k) {
            int sk = spinBasis[k % subDim];
            int hk = holeBasis[int(k/subDim)];
            for (int l = 0; l < len; ++l) {
                int sl = spinBasis[l % subDim];
                int hl = holeBasis[int(l/subDim)];
                if (j == hk && i == hl && sk == sl) {
                    res += pow(-1, hk+hl+1)*std::conj(v[k])*v[l];
                    }
                }
            }
        return res;
        }

// Change the groud state into the Marshall sign basis.
template<typename T>
void tjChainHalf<T>::Marshall(T* v1, T* v2) { 
        int len = Dim();
        for (int i = 0; i < len; ++i) {
            int s = spinBasis[i % subDim];
            int h = holeBasis[int(i/subDim)];
            std::bitset<numSite> spinConfig(s);
            int n = 0;
            for (int j = 0; j < numSite; ++j){
                if ( 0 == j % 2) {
                    if (j == h) {++n;} // Only for the electron removed with spin down.
                    else if (j < h && 0 == spinConfig[j]) { ++n; }
                    else if (j > h && 0 == spinConfig[j-1]) { ++n; }
                    }
                }
            v2[i] = pow(-1, n)*pow(-1, h)*v1[i]; // The second pow() only exists for the electron removed with spin down.
            }
        }

