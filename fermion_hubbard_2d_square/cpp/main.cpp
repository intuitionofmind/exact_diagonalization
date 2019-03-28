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
        int numSiteX = 2;
        int numSiteY = 4;
        double U = 0.0;
        int numSam = 100;
        int numEval = 10;

        std::ofstream file_data("results.dat", std::ios_base::app);
        file_data << "u," << "gsE," <<  "trans" << std::endl;

        for (int s = 0; s < numSam; ++s) {
            Hubbard<std::complex<double>> Hub(U, numEleUp, numEleDown, numSiteX, numSiteY, "OBC", "OBC");
            int dim = Hub.HilbertDim();
            // std::cout << dim << std::endl;
            if (0 == s) { Hub.SaveHilbert(); }
            ARCompStdEig<double, Hubbard<std::complex<double>>> EigProb;
            EigProb.DefineParameters(dim, numEval, &Hub, &Hubbard<std::complex<double>>::Hamiltonian, "SR", 30);
            auto eigVal = new std::complex<double>[numEval];
            auto eigVec = new std::complex<double>[numEval*dim];
            int nconv = EigProb.EigenValVectors(eigVec, eigVal);

            // To sort the eigenvalues and keep track of the index.
            std::vector<std::pair<double, int>> temp;
            for (int j = 0; j < nconv; ++j) { temp.push_back(std::pair<double, int> (std::real(eigVal[j]), j)); }
            std::sort(temp.begin(), temp.end());
            // for (auto a : temp) { std::cout << a.first << " " << a.second << std::endl; }
            // std::cout << "eval: " << std::setprecision(14) << temp[j] / (numSiteX*numSiteY) << std::endl;

            for (int j = 0; j < nconv; ++j) {
                double val = temp[j].first;
                int index = temp[j].second;
                file_eigvals.write((char*)(&val), sizeof(double));

                for (int k = 0; k < dim; ++k) {
                    int d = index*dim+k;
                    double rr = std::real(eigVec[d]);
                    double ii = std::imag(eigVec[d]);
                    file_eigvecs.write((char*)(&rr), sizeof(double));
                    file_eigvecs.write((char*)(&ii), sizeof(double));
                    }
                }

            auto gsWf = new std::complex<double>[dim];
            auto gsWfTemp = new std::complex<double>[dim];
            auto gsWfTempTemp = new std::complex<double>[dim];
            for (int k = 0; k < dim; ++k) {
                gsWf[k] = eigVec[(temp[0].second)*dim+k];
                }
            Hub.TotalSPlusSMinus(gsWf, gsWfTemp);
            Hub.TotalSMinusSPlus(gsWf, gsWfTempTemp);
            auto transv = 0.5*(Hub.Dot(gsWf, gsWfTemp)+Hub.Dot(gsWf, gsWfTempTemp));

            delete [] eigVal;
            delete [] eigVec;
            delete [] gsWf;
            delete [] gsWfTemp;
            delete [] gsWfTempTemp;

            std::cout << s << " " <<std::setprecision(4) << U << " " << std::setprecision(8) << temp[0].first << " " << temp[1].first << " " << transv << std::endl;
            file_data << std::setprecision(4) << U << "," << std::setprecision(8) << temp[0].first << "," << std::real(transv) << std::endl;
            U += 2.0;

            }
        end = time(NULL);
        file_log << "Time: " << (end-start)/60.0 << " min" << std::endl;
        std::cout << "Time: " << (end-start) << " s" << std::endl;
file_eigvals.close();
        file_eigvecs.close();
        file_log.close();
        return 1;
    }
