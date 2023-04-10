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

    // Neutron PDFs (using isospin symmetry)
    double u_n = d_p;
    double d_n = u_p;

    // Deuteron PDFs (using isospin symmetry)
    double u_d = 0.5 * (u_p + d_n);
    double d_d = 0.5 * (d_p + u_n);

    // PDF ratios
    double F_pu = u_p / (u_p + d_p);
    double F_pd = d_p / (u_p + d_p);
    double F_du = u_d / (u_d + d_d);
    double F_dd = d_d / (u_d + d_d);

    std::cout << "F_pu: " << F_pu << std::endl;
    std::cout << "F_pd: " << F_pd << std::endl;
    std::cout << "F_du: " << F_du << std::endl;
    std::cout << "F_dd: " << F_dd << std::endl;

    delete pdf;
    return 0;
}
