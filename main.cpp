#include <iostream>
#include <vector>

class Multigrid {
    int depth, iter, N;
    double x0, x1, u0, u1;
    bool verbose;
    std::vector<int> num_cells;
    std::vector<double> dx;
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> s;
    std::vector<std::vector<double>> res;
    Multigrid(int depth, int iter, double x0, double x1, double u0, double u1, std::vector<double> s, bool verbose);
    void relax(int depth);
    void relax_rb(int depth);
    void prolongation(int depth);
    void restriction(int depth);
    void residual(int depth);
    void multigrid(int lvs);
};

Multigrid::Multigrid(int depth, int iter, double x0, double x1, double u0, double u1, std::vector<double> s, bool verbose) {
    this -> depth = depth;
    this -> N = 1 << depth;
    this -> iter = iter;
    this -> x0 = x0;
    this -> x1 = x1;
    this -> u0 = u0;
    this -> u1 = u1;
    this -> verbose = verbose;
    for (int i = 0; i < depth; ++i) {
        this -> dx.emplace_back((x1 - x0) / (2 << i));
        this -> num_cells.emplace_back(2 << i);
        std::vector<double> now((2 << i) + 1);
        this -> x.emplace_back(now);
        this -> u.emplace_back(now);
        this -> s.emplace_back(now);
        this -> res.emplace_back(now);
        for (int j = 0; j < ((2 << i) + 1); ++j) {
            this -> x[i][j] = x0 + j * this -> dx[i];
        }
    }

    this -> u[depth - 1][0] = u0;
    this -> u[depth - 1][1 << depth] = u1;
    for (int i = 1; i < this -> N; ++i) {
        this -> s[depth - 1][i] = s[i];
    }
};

int main() {
    std::vector<std::vector<int>> t;
    std::vector<int> a(3, 1);
    t.emplace_back(a);
    t.emplace_back(a);
    t[1][0] = 2;
    std::cout << t[0][0] << " " << t[1][0] << "\n";
    return 0;
}
