#include <iostream>
#include <cmath>
#include <LHAPDF/LHAPDF.h>

int main() {
    double x = 0.1;
    double Q2 = 10.0;
    double Q = std::sqrt(Q2);

    LHAPDF::PDF* pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);

    double u_p = pdf->xfxQ(2, x, Q);
    double d_p = pdf->xfxQ(1, x, Q);

    std::cout << "u_p: " << u_p << std::endl;
    std::cout << "d_p: " << d_p << std::endl;

    delete pdf;
    return 0;
}
