#include "head.hpp"
#include "class_hubbard.hpp"

int main() {
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);

        time_t start, end;
        start = time(NULL);

        int numEleUp = 4;
        int numEleDown = 4;
        int numSiteX = 8;
        int numSiteY = 1;
        double U = 1.0;
        int numSam = 1;
        int numEval = 10;

        // PrintHam(J);
        // file_log << numSite << std::endl;
        // file_log << numSam << std::endl;
        // file_log << numEval << std::endl;
        // file_log << dim << std::endl;
        // file_log << U << std::endl;
        // file_log << step << std::endl;

        // Test.
        // int n = 10;;
        // std::bitset<4> bb(7);
        // boost::dynamic_bitset<> bb(n, 2);
        // for (int j = 0; j < n; ++j) { std::cout << bb[j]; }
        // std::cout << std::endl;

        Hubbard<arcomplex<double>> Hub(U, numEleUp, numEleDown, numSiteX, numSiteY, "PBC", "OBC");
        int dim = Hub.HilbertDim();
        auto u = new arcomplex<double>[dim];
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        Hub.SetOne(u, 2);
        Hub.SetOne(w, 0);
        Hub.Hamiltonian(u, v);
//        std::cout << "Matrix element: "  << Hub.Dot(u, v) << std::endl;
        delete [] u;
        delete [] v;
        delete [] w;
        // for (int j = 0; j < dim; ++j) { v[j] = 0.1; }
        // Hub.Hamiltonian(v, w);
        // Hub.PrintHam();


        for (int i = 0; i < numSam; ++i) {
            // std::cout << i << std::endl;
            ARCompStdEig<double, Hubbard<arcomplex<double>>> prob;
            prob.DefineParameters(dim, numEval, &Hub, &Hubbard<arcomplex<double>>::Hamiltonian, "SR");
        
            auto eigVal = new arcomplex<double>[numEval];
            auto eigVec = new arcomplex<double>[numEval*dim];
            int nconv = prob.EigenValVectors(eigVec, eigVal);
            // std::vector<int> order;
            // H.SortEval(nconv, eigVal, eigValSort, order);
            for (int j = 0; j < nconv; ++j) {
                double r = std::real(eigVal[j]);
                file_eigvals.write((char*)(&r), sizeof(double));
                // if (0 == i) { std::cout << "eval: " << std::setprecision(14) << eigVal[j] << std::endl; }
                if (0 == i) { std::cout << "eval: " << std::setprecision(14) << std::real(eigVal[j]) / (numSiteX*numSiteY) << std::endl; }
                // for (int k = 0; k < dim; ++k) {
                    // v1[k] = eigVec[order[j]*dim+k];
                    // }
                for (int k = 0; k < dim; ++k) {
                    double rr = std::real(eigVec[j*dim+k]);
                    double ii = std::imag(eigVec[j*dim+k]);
                   // if (0 == i && 0 == j) { std::cout << rr << " " << ii << std::endl; }
                    file_eigvecs.write((char*)(&rr), sizeof(double));
                    file_eigvecs.write((char*)(&ii), sizeof(double));
                    }
                }
            delete [] eigVal;
            delete [] eigVec;
            file_log << (i+1) << std::endl;
            U += 0.1;
            }

        end = time(NULL);
        file_log << "Time: " << (end-start)/60.0 << " min" << std::endl;
        std::cout << "Time: " << (end-start) << " s" << std::endl;
        file_eigvals.close();
        file_eigvecs.close();
        file_log.close();
        return 1;
        }
