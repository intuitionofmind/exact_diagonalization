template<class T>
class heisenHalf {
        private:
        int n;

        public:
        int Dim();
        void MultVec(T* v, T* w);
        T Dot(T* v1, T* v2);
        void SetOne(T* v, int i);
        heisenHalf(int d) { n = d; }  // constructor
        };

template<class T>
int heisenHalf<T>::Dim() { return n; }

/*
One dimensional spin-1/2 Heisenberg Hamiltonian including a 1/4 subtraction.

Required by Arpackpp package:
There only requirements make by ARPACK++ are that MultVec musth have two pointers to vectors of type T as paraments and the input vector must precede the output vector.
*/
template<class T>
void heisenHalf<T>::MultVec(T* v, T* w) {
        int len = Dim();
        for (int i = 0; i < len; ++i) { w[i] = 0.0; }
        for (int i = 0; i < len; ++i) {
            int s = bas[i];
            std::bitset<numSite> config(s);
            // std::cout << i << << config << std::endl;
            int cyc = 0;
            if ("PBC" == flagBoun) { cyc = numSite; }
            if ("OBC" == flagBoun) { cyc = numSite-1; }
            for (int j = 0; j < cyc; ++j) {
                int jj = (j+1) % numSite;
                if (config[j] == config[jj]) { continue; }
                else {
                    w[i] -= 0.5*v[i];
                    std::bitset<numSite> temp (config);
                    temp.flip(j);
                    temp.flip(jj);
                    int ss = temp.to_ulong();
                    std::vector<int>::iterator k = std::lower_bound(bas.begin(), bas.end(), ss);
                    w[std::distance(bas.begin(), k)] += 0.5*v[i];
                    }
                }
            }
        }

template<typename T>
T heisenHalf<T>::Dot(T* v1, T* v2) {
        int len = Dim();
        T res = 0.;
        for (int i = 0; i < len; i++) { res += v1[i]*v2[i]; }
        return res;
        }

template<typename T>
void heisenHalf<T>::SetOne(T* v, int i) {
        int len = Dim();
        for (int j = 0; j < len; j++) { v[j] = 0.; }
        v[i] = 1.;
        }
/*
int Flip(int s, int i, int j) {
        int f = int(pow(2, i))+int(pow(2, j));
        return (s ^ f);
        }
*/
