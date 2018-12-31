#include "head.hpp"
#include "class_fermion.hpp"

int main() {
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);

        time_t start, end;
        start = time(NULL);

        int numEle = 4;
        int numSiteX = 8;
        int numSiteY = 1;
        int numEval = 20;

        /*
        FreeFermion<arcomplex<double>> F(numEle, numSiteX, numSiteY, "PBC", "OBC");
        int dim = F.HilbertDim();
        auto u = new arcomplex<double>[dim];
        auto v = new arcomplex<double>[dim];
        auto w = new arcomplex<double>[dim];
        F.SetOne(u, 2);
        F.SetOne(w, 0);
        F.Hamiltonian(u, v);
        std::cout << "Matrix element: "  << F.Dot(u, v) << std::endl;
        delete [] u;
        delete [] v;
        delete [] w;
        // F.PrintHam();
        */

        FreeFermion<arcomplex<double>> F(numEle, numSiteX, numSiteY, "PBC", "OBC");
        int dim = F.HilbertDim();
        // std::cout << i << std::endl;
        ARCompStdEig<double, FreeFermion<arcomplex<double>>> prob;
        prob.DefineParameters(dim, numEval, &F, &FreeFermion<arcomplex<double>>::Hamiltonian, "SR");

        auto eigVal = new arcomplex<double>[numEval];
        auto eigVec = new arcomplex<double>[numEval*dim];
        int nconv = prob.EigenValVectors(eigVec, eigVal);
        for (int j = 0; j < nconv; ++j) {
            double r = std::real(eigVal[j]);
            file_eigvals.write((char*)(&r), sizeof(double));
            std::cout << "eval: " << std::setprecision(14) << std::real(eigVal[j]) / (numSiteX*numSiteY) << std::endl;
            for (int k = 0; k < dim; ++k) {
                double rr = std::real(eigVec[j*dim+k]);
                double ii = std::imag(eigVec[j*dim+k]);
                file_eigvecs.write((char*)(&rr), sizeof(double));
                file_eigvecs.write((char*)(&ii), sizeof(double));
                }
            }
        delete [] eigVal;
        delete [] eigVec;

        end = time(NULL);
        std::cout << "Time: " << (end-start) << std::endl;
        file_eigvals.close();
        file_eigvecs.close();
        file_log.close();
        return 1;
    }
