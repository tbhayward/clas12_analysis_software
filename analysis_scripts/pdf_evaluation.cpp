#include <iostream>
#include <cmath>
#include <vector>
#include <LHAPDF/LHAPDF.h>

int main() {
    // ppip
    std::vector<double> meanXRGA = {0.294096, 0.282596, 0.260506, 0.238092, 0.215682, 0.192768};
    std::vector<double> meanQ2RGA = {2.99703, 2.88093, 2.72367, 2.54051, 2.37417, 2.22758};
    std::vector<double> meanXRGB = {0.281038,0.274911,0.258215,0.24149,0.224762,0.208944};
    std::vector<double> meanQ2RGB = {2.98496,2.9152,2.79621,2.66012,2.52836,2.43014};

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

    std::vector<double> R_u_p_values, R_d_p_values, R_u_d_values, R_d_d_values;

    for (size_t i = 0; i < meanXRGA.size(); i++) {

        double u_p = pdf->xfxQ(2, meanXRGA[i], std::sqrt(meanQ2RGA[i]));
        double d_p = pdf->xfxQ(1, meanXRGA[i], std::sqrt(meanQ2RGA[i]));

        double R_u_p = u_p / (u_p + d_p);
        double R_d_p = d_p / (u_p + d_p);

        double u_n = pdf->xfxQ(-1, meanXRGB[i], std::sqrt(meanQ2RGB[i]));
        double d_n = pdf->xfxQ(-2, meanXRGB[i], std::sqrt(meanQ2RGB[i]));
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
