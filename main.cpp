#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define PI 3.141592653589793238462

class Multigrid {
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
    Multigrid(int depth, int iter, double x0, double x1, double u0, double u1, std::vector<double> s, bool verbose);
    void relax(int depth);
    void relax_rb(int depth);
    void prolongation(int depth);
    void restriction(int depth);
    void residual(int depth);
    void multigrid(int lvs);
    std::vector<double> solve();
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

void Multigrid::relax(int depth) {
    for (int i = 0; i < this -> iter; ++i) {
        this -> relax_rb(depth);
    }
}

void Multigrid::relax_rb(int depth) {
    for (int i = 1; i < this -> num_cells[depth]; i+=2) {
        this -> u[depth][i] = (this -> u[depth][i + 1] + this -> u[depth][i - 1] - this -> s[depth][i] * this -> dx[depth] * this -> dx[depth]) / 2;
    }
    for (int i = 2; i < this -> num_cells[depth]; i+=2) {
        this -> u[depth][i] = (this -> u[depth][i + 1] + this -> u[depth][i - 1] - this -> s[depth][i] * this -> dx[depth] * this -> dx[depth]) / 2;
    }
}

void Multigrid::prolongation(int depth) {
    for (int i = 0; i < this -> num_cells[depth]; ++i) {
        if (i & 1) {
            this -> u[depth][i] += (this -> u[depth - 1][i >> 1] + this -> u[depth - 1][(i >> 1) + 1]) / 2;
        }
        else {
            this -> u[depth][i] += this -> u[depth - 1][i >> 1];
        }
    }
}

void Multigrid::restriction(int depth) {
    this -> s[depth][0] = this -> res[depth + 1][0];
    this -> s[depth][s[depth].size() - 1] = this -> res[depth + 1][res[depth + 1].size() - 1];
    for (int i = 1; i < this -> num_cells[depth]; ++i) {
        this -> s[depth][i] = (this -> res[depth + 1][2 * i - 1] + 2 * this -> res[depth + 1][2 * i] + this -> res[depth + 1][2 * i + 1]) / 4;
    }
    for (int i = 0; i < this -> num_cells[depth]; ++i) {
        this -> u[depth][i] = 0;
    }
}

void Multigrid::residual(int depth) {
    for (int i = 1; i < this -> num_cells[depth]; ++i) {
        this -> res[depth][i] = this -> s[depth][i] - (this -> u[depth][i + 1] + this -> u[depth][i - 1] - 2 * this -> u[depth][i]) / this -> dx[depth] / this -> dx[depth];
    }
}

void Multigrid::multigrid(int lvs) {
    if (lvs == 0) {
        this -> relax(lvs);
        return;
    }
    this -> relax(lvs);
    this -> residual(lvs);
    this -> restriction(lvs - 1);
    this ->multigrid(lvs - 1);
    this ->prolongation(lvs);
    this ->relax(lvs);
}

std::vector<double> Multigrid::solve() {
    for (int i = 0; i < 10; ++i) {
        this ->multigrid(this -> depth - 1);
    }

    return this -> u[this -> u.size() - 1];
}

double src_func(double x) {
    return - PI * PI * std::sin(PI * x);
}

int main() {
    int level = 19;
    double x0 = 0.0;
    double x1 = 1.0;
    double u0 = 1.0;
    double u1 = 0.0;
    int N = 1 << level;
    double dx = (x1 - x0) / N;
    std::vector<double> x(N + 1);
    std::vector<double> source(N + 1);
    for (int i = 0; i < N + 1; ++i) {
        x[i] = x0 + i * dx;
    }
    for (int i = 0; i < N + 1; ++i) {
        source[i] = src_func(x[i]);
    }

    Multigrid mg(level, 1, x0, x1, u0, u1, source, false);

    std::vector<double> sol = mg.solve();

    std::ofstream fout;
    fout.open("sol.txt");
    if (!fout) return 0;
    for (double i : x) {
        fout << i << " ";
    }
    fout << "\n";
    for (double i : sol) {
        fout << i << " ";
    }

    fout.close();
    return 0;
}
