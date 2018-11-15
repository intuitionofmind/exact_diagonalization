#include "head.h"
#include "precondition.hpp"
#include "hamiltonian.hpp"
#include "arssym.h"

int main() {
        std::ofstream file_log("log", std::ios_base::app);
        std::ofstream file_eigvals("eigenvalues.dat", std::ios_base::app | std::ios_base::binary);
        std::ofstream file_eigvecs("eigenvectors.dat", std::ios_base::app | std::ios_base::binary);
        auto eigVal = new double[numEval];
        auto eigVec = new double[(numEval+1)*dim];
 
        time_t start, end;
        // PrintHam();
        start = time(NULL);
        heisenHalf<double> H(dim);
        ARSymStdEig<double, heisenHalf<double>> prob;
        prob.DefineParameters(dim, numEval, &H, &heisenHalf<double>::MultVec, "SA");
        int nconv = prob.EigenValVectors(eigVec, eigVal);
        for (int i = 0; i < nconv; i++) {
            std::cout << (eigVal[i]/(numSite)+0.25)*numSite << std::endl; // The exact energy is -0.4438.
            // std::cout << eigVal[i] << std::endl; // The exact energy is -0.4438.
            for (int j = 0; j < dim; j++) {
                file_eigvecs.write((char*)(&eigVec[i*dim+j]), sizeof(double));
                }
            }

        delete [] eigVal;
        delete [] eigVec;
        end = time(NULL);
        file_log << "Time: " << double(end-start)/60.0 << " min" << std::endl;
        file_log.close();
        file_eigvals.close();
        file_eigvecs.close();
        return 1;
        }

/* 
PrintHam() function is used to print the Hamiltonian matrix explicitly to check whether the code is right or not.
*/
void PrintHam() {
        heisenHalf<double> A(dim);
        double* v1 = new double[dim];
        double* v2 = new double[dim];
        double* temp = new double[dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                A.SetOne(v1, i);
                A.SetOne(v2, j);
                A.MultVec(v2, temp);
                std::cout << A.Dot(v1, temp) << " ";
                }
            std::cout << std::endl;
            }
        delete [] v1;
        delete [] v2;
        delete [] temp;
        }
