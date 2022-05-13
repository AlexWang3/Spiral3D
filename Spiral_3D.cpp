#include <iostream>
#include <fstream>
#include <vector>
#include <math.h> 
using namespace std;
vector<vector<double>> V;
vector<vector<int>> F;

const double e = 2.71828;
double a = 1;
double b = log(20) * -1;

int output_file() {
    ofstream File("spiral.obj");
    for (int i = 0; i < V.size(); i++) {
        File << "v";
        for (int j = 0; j < 3; j++) {
            File << " ";
            File << V[i][j];
        }
        File << "\n";
    }
    for (int i = 0; i < F.size(); i++) {
        File << "f";
        for (int j = 0; j < 3; j++) {
            File << " ";
            File << F[i][j];
        }
        File << "\n";
    }
    File.close();
    return 0;
}

vector<double> cexp(vector<double> p) {
    double r = pow(e, p[0]);
    return { r * cos(p[1]), r * sin(p[1]) };
}

vector<double> cadd(vector<double> p, vector<double> q) {
    return { p[0] + q[0], p[1] + q[1] };
}

vector<double> cmul(vector<double> p, vector<double> q) {
    return { p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0] };
}

vector<double> crecip(vector<double> p)
{
    double d = p[0] * p[0] + p[1] * p[1];
    return { p[0] / d, -p[1] / d };
}

vector<double> cdiv(vector<double> p, vector<double> q) {
    return cmul(p, crecip(q));
}

vector<double> m(double t, double dy) {
    vector<double> P = { t * a, t * b + dy };
    vector<double> ep = cexp(P);
    return cdiv(cadd(ep, { -1, 0 }), cadd(ep, { 1, 0 }));
}

double map(int idx, int start, int end, double target_start, double target_end) {
    double part = ((double)idx - (double)start) / ((double)end - (double)start);
    return target_start + part * (target_end - target_start);
}

double b_start = 20.0;
double b_end = 100.0;
double q_min = -1.0;
double q_max = 1.0;
int layer_points = 5000;
double layer_height = 2.0;
int layer_num = 50;
int layer_spiral_num = 2;

void drawThing(vector<vector<vector<double>>>& spirals, double qr) {

    vector<vector<double>> s1;
    vector<vector<double>> s2;

    for (int idx = 0; idx < layer_points; idx++) {
        double t = map(idx, 0, layer_points, -100, 100);
        vector<double> P = m(t, qr);
        s1.push_back({ P[0] * 40, P[1] * 40 });
    }

    for (int idx = 0; idx < layer_points; idx++) {
        double t = map(idx, 0, layer_points, -100, 100);
        vector<double> P = m(t, -qr);
        s2.push_back({ P[0] * 40, P[1] * 40 });
    }
    spirals.push_back(s1);
    spirals.push_back(s2);
}

void update_a_b(double b_val) {
    a = 1;
    b = log(b_val) * -1;
    double M = max(a, b);
    a /= M;
    b /= M;
}

int main()
{

    // Create pair spirals
    double b_val = b_start;
    update_a_b(b_val);
    double b_diff = (b_end - b_start) / (layer_num - 1.0);
    double qr = q_min;
    double q_diff = (q_max - q_min) / (layer_num - 1.0);

    /*
    vector<vector<double>> intialSpiral;
    for (int idx = 0; idx < layer_points; idx++) {
        double t = map(idx, 0, layer_points, -100, 100);
        vector<double> P = m(t, qr);
        intialSpiral.push_back({ P[0] * 40, P[1] * 40 });
    }
    */

    vector<vector<vector<double>>> spirals;
    for (int i = 0; i < layer_num; i++) {
        drawThing(spirals, qr);
        qr += q_diff;
        b_val += b_diff;
        update_a_b(b_val);
    }
    
    /*
    for (int i = 0; i < layer_points; i++) {
        V.push_back({ intialSpiral[i][0], 0, intialSpiral[i][1] });
    }
    */

    double cur_height = 0;
    for (int i = 0; i < layer_num; i++) {
        for (int j = 0; j < layer_spiral_num; j++) {
            vector<vector<double>> cur_spiral = spirals[i * layer_spiral_num + j];
            for (int p = 0; p < layer_points; p++) {
                V.push_back({ cur_spiral[p][0], cur_height, cur_spiral[p][1] });
            }
        }
        cur_height -= layer_height;
    }

    
    for (int i = 0; i < spirals.size() - layer_spiral_num; i++) {
        for (int j = 1 + i * layer_points; j < (i + 1) * layer_points; j++) {
            F.push_back({ j, j + 1, j + layer_spiral_num * layer_points });
            F.push_back({ j + 1, j + layer_spiral_num * layer_points, j + 1 + layer_spiral_num * layer_points });
        }
    }

    for (int i = 1; i < layer_points; i++) {
        F.push_back({ i, i + 1, i + layer_points });
        F.push_back({ i + 1, i + layer_points, i + 1 + layer_points });
    }

    int start_of_last_group = (spirals.size() - layer_spiral_num) * layer_points;

    for (int i = 1 + start_of_last_group; i < start_of_last_group + layer_points; i++) {
        F.push_back({ i, i + 1, i + layer_points });
        F.push_back({ i + 1, i + layer_points, i + 1 + layer_points });
    }

    output_file();
}
