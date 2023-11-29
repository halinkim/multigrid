//
// Created by halin on 2023-11-29.
//

#ifndef MULTIGRID_MULTIGRID3D_H
#define MULTIGRID_MULTIGRID3D_H


#include <vector>

class Multigrid3d {
public:
    int depth, iter, N;
    double x0, x1, u0, u1;
    bool verbose;
    std::vector<int> num_cells;
    std::vector<double> dx;
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> s;
    std::vector<std::vector<double>> res;
    Multigrid3d(int depth, int iter, double x0, double x1, double u0, double u1, std::vector<double> s, bool verbose);
    void relax(int depth);
    void relax_rb(int depth);
    void prolongation(int depth);
    void restriction(int depth);
    void residual(int depth);
    void multigrid(int lvs);
    std::vector<double> solve();
};


#endif //MULTIGRID_MULTIGRID3D_H
