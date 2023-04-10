#include <iostream>
#include <cmath>
#include <vector>
#include <LHAPDF/LHAPDF.h>

int main() {
    std::vector<double> meanX = {0.294096, 0.282596, 0.260506, 0.238092, 0.215682, 0.192768};
    std::vector<double> meanQ2 = {2.99703, 2.88093, 2.72367, 2.54051, 2.37417, 2.22758};

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

    std::vector<double> R_u_p_values, R_d_p_values, R_u_d_values, R_d_d_values;

    for (size_t i = 0; i < meanX.size(); i++) {
        double x = meanX[i];
        double Q2 = meanQ2[i];

        double u_p = pdf->xfxQ(2, x, std::sqrt(Q2));
        double d_p = pdf->xfxQ(1, x, std::sqrt(Q2));

        double R_u_p = u_p / (u_p + d_p);
        double R_d_p = d_p / (u_p + d_p);

        double u_n = pdf->xfxQ(-1, x, std::sqrt(Q2));
        double d_n = pdf->xfxQ(-2, x, std::sqrt(Q2));
        double R_u_d = 0.5 * (u_p + u_n) / (u_p + d_p + u_n + d_n);
        double R_d_d = 0.5 * (d_p + d_n) / (u_p + d_p + u_n + d_n);

        R_u_p_values.push_back(R_u_p);
        R_d_p_values.push_back(R_d_p);
        R_u_d_values.push_back(R_u_d);
        R_d_d_values.push_back(R_d_d);
    }

    delete pdf;

    auto print_vector = [](const std::vector<double>& v, const std::string& name) {
        std::cout << name << " = {";
        for (size_t i = 0; i < v.size(); i++) {
            std::cout << v[i];
            if (i < v.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "}" << std::endl;
    };

    print_vector(R_u_p_values, "R_u_p");
    print_vector(R_d_p_values, "R_d_p");
    print_vector(R_u_d_values, "R_u_d");
    print_vector(R_d_d_values, "R_d_d");

    return 0;
}
