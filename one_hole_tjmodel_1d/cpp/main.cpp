#include "head.h"
#include "precondition.hpp"
#include "hamiltonian.hpp"
#include "arssym.h"
#include "arcomp.h"
#include "arscomp.h"

int main() {
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);
        time_t start, end;
        start = time(NULL);
        auto eigVal = new arcomplex<double>[numEval];
        auto eigValSort = new arcomplex<double>[numEval];
        auto eigVec = new arcomplex<double>[numEval*dim];
        auto v1 = new arcomplex<double>[dim];
        auto v2 = new arcomplex<double>[dim];
        double J = 0.3;
        PrintHam(40.3);
        file_log << numSite << std::endl;
        file_log << numSam << std::endl;
        file_log << numEval << std::endl;
        file_log << dim << std::endl;
        file_log << J << std::endl;
        file_log << step << std::endl;
        file_log << sigma << std::endl;

        for (int i = 0; i < numSam; ++i) {
            std::cout << i << std::endl;
            tjChainHalf<arcomplex<double>> H(dim, J);
            ARCompStdEig<double, tjChainHalf<arcomplex<double>>> prob;
            prob.DefineParameters(dim, numEval, &H, &tjChainHalf<arcomplex<double>>::MultVec, "SR");
            int nconv = prob.EigenValVectors(eigVec, eigVal);
            std::vector<int> order;
            H.SortEval(nconv, eigVal, eigValSort, order);
            for (int j = 0; j < nconv; ++j) {
                double r = std::real(eigValSort[j]);
                file_eigvals.write((char*)(&r), sizeof(double));
                // if (0 == i) { std::cout << "eval: " << std::setprecision(14) << (eigValSort[j])*3.0 << std::endl; }
                if (0 == i) { std::cout << "eval: " << std::setprecision(14) << (eigValSort[j]) << std::endl; }
                for (int k = 0; k < dim; ++k) {
                    v1[k] = eigVec[order[j]*dim+k];
                    }
                for (int k = 0; k < dim; ++k) { v2[k] = v1[k]; }
                for (int k = 0; k < dim; ++k) {
                    double rr = std::real(v1[k]);
                    double ii = std::imag(v1[k]);
                    // if (0 == i && 0 == j) { std::cout << rr << " " << ii << std::endl; }
                    file_eigvecs.write((char*)(&rr), sizeof(double));
                    file_eigvecs.write((char*)(&ii), sizeof(double));
                    }
                }
            file_log << (i+1) << std::endl;
            J += 0.1;
            }
        delete [] eigVal;
        delete [] eigValSort;
        delete [] eigVec;
        end = time(NULL);
        file_log << "Running time: " << (end-start)/60.0 << " min" << std::endl;
        file_log.close();
        return 1;
    }

        // PrintHam() function is used to print the Hamiltonian matrix explicitly to check whether the code is right or not.
void PrintHam(double JJ) {
        tjChainHalf<double> H(dim, JJ);
        std::ofstream file_matrix("matrix.txt", std::ios_base::app);
        auto v1 = new double[dim];
        auto v2 = new double[dim];
        auto temp = new double[dim];
        std::cout << dim << std::endl;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                for (int l = 0; l < dim; ++l) {
                    v1[l] = 0.0;
                    v2[l] = 0.0;
                    temp[l] = 0.0;
                    }
                v1[i] = 1.0;
                v2[j] = 1.0;
                H.MultVec(v2, temp);
                file_matrix << H.Dot(v1, temp) << std::endl;
                // std::cout << i*dim+j << std::endl;
                }
            }
        delete [] v1;
        delete [] v2;
        delete [] temp;
        file_matrix.close();
        }
