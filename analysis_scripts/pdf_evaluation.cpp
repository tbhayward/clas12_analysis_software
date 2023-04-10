#include <iostream>
#include <cmath>
#include <LHAPDF/LHAPDF.h>

int main(int argc, char* argv[]) {
    std::cout << std::endl << std::endl;
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " x Q2" << std::endl;
        return 1;
    }

    double x = std::stod(argv[1]);
    double Q2 = std::stod(argv[2]);

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

    double u_p = pdf->xfxQ(2, x, std::sqrt(Q2));
    double d_p = pdf->xfxQ(1, x, std::sqrt(Q2));

    double R_u_p = u_p / (u_p + d_p);
    double R_d_p = d_p / (u_p + d_p);

    double u_n = pdf->xfxQ(2, x, std::sqrt(Q2), 1);
    double d_n = pdf->xfxQ(1, x, std::sqrt(Q2), 1);
    double R_u_d = 0.5 * (u_p + u_n) / (u_p + d_p + u_n + d_n);
    double R_d_d = 0.5 * (d_p + d_n) / (u_p + d_p + u_n + d_n);


    std::cout << "R_u_p: " << R_u_p << std::endl;
    std::cout << "R_d_p: " << R_d_p << std::endl;
    std::cout << "R_u_d: " << R_u_d << std::endl;
    std::cout << "R_d_d: " << R_d_d << std::endl;

    delete pdf;
    return 0;
}
