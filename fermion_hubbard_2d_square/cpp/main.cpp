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
        double U = 20.0;
        // int numSam = 1;
        int numEval = 4;

        Hubbard<std::complex<double>> Hub(U, numEleUp, numEleDown, numSiteX, numSiteY, "PBC", "OBC");
        int dim = Hub.HilbertDim();
        // Hub.PrintHam();

        ARCompStdEig<double, Hubbard<std::complex<double>>> prob;
        prob.DefineParameters(dim, numEval, &Hub, &Hubbard<std::complex<double>>::Hamiltonian, "SR");
 
        auto eigVal = new std::complex<double>[numEval];
        auto eigVec = new std::complex<double>[numEval*dim];
        int nconv = prob.EigenValVectors(eigVec, eigVal);

        std::vector<double> temp;
        for (int j = 0; j < nconv; ++j) { temp.push_back(std::real(eigVal[j])); }
        std::sort(temp.begin(), temp.end());

        for (int j = 0; j < nconv; ++j) {
            double r = std::real(eigVal[j]);
            file_eigvals.write((char*)(&r), sizeof(double));
            // std::cout << "eval: " << std::setprecision(14) << temp[j] / (numSiteX*numSiteY) << std::endl;
            std::cout << "eval: " << std::setprecision(14) << temp[j] << std::endl;
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
        file_log << "Time: " << (end-start)/60.0 << " min" << std::endl;
        std::cout << "Time: " << (end-start) << " s" << std::endl;
        file_eigvals.close();
        file_eigvecs.close();
        file_log.close();
        return 1;
        }
