#include <iostream>
#include <cmath>
#include <LHAPDF/LHAPDF.h>

int main(int argc, char* argv[]) {
    std::cout << endl << endl;
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " x Q2" << std::endl;
        return 1;
    }

    double x = std::stod(argv[1]);
    double Q2 = std::stod(argv[2]);
    double Q = std::sqrt(Q2);

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

    double u_p = pdf->xfxQ(2, x, Q);
    double d_p = pdf->xfxQ(1, x, Q);

    std::cout << "u_p: " << u_p << std::endl;
    std::cout << "d_p: " << d_p << std::endl;
    std::cout << endl << endl;

    delete pdf;
    return 0;
}
