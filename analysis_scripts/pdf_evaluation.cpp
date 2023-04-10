#include <iostream>
#include <cmath>
#include <vector>
#include <LHAPDF/LHAPDF.h>

int main() {
    std::vector<double> meanX = {0.294096, 0.282596, 0.260506, 0.238092, 0.215682, 0.192768};
    std::vector<double> meanQ2 = {2.99703, 2.88093, 2.72367, 2.54051, 2.37417, 2.22758};

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

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

        std::cout << "x = " << x << ", Q2 = " << Q2 << std::endl;
        std::cout << "R_u_p: " << R_u_p << std::endl;
        std::cout << "R_d_p: " << R_d_p << std::endl;
        std::cout << "R_u_d: " << R_u_d << std::endl;
        std::cout << "R_d_d: " << R_d_d << std::endl;
        std::cout << std::endl;
    }

    delete pdf;
    return 0;
}
