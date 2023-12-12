
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

#include"src/pr/pr_imp.h"

namespace thermolib {

void PRimp::SetFluid(double Tc, double pc, double R, double omega) {
    Tc_ = Tc;
    pc_ = pc;
    R_ = R;
    kappa_ = 0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2);
    ac_ = 0.45724 * pow(R_, 2) * pow(Tc_, 2) / pc_;
    bc_ = 0.07780 * R_ * Tc_ / pc_;
    Zc_ = 0.307;
}

void PRimp::TPflash(double T, double p) {
    is_single_phase_ = true;  // TPflash 得到单相区
    T_ = T;
    P_ = p;
    double a = ac_ * pow(1 + kappa_ * (1 - sqrt(T_ / Tc_)), 2);
    double b = bc_;
    A_ = a * P_ / pow(R_, 2) / pow(T_, 2);
    B_ = b * P_ / R_ / T_;
    double Zg = 0, Zl = 0;  // 用于记录一元三次方程求解结果
    if (ShengjinFanSolver(1, -(1 - B_), A_ - 3 * pow(B_, 2) - 2 * B_, -(A_ * B_ - pow(B_, 2) - pow(B_, 3)), &Zg, &Zl) == true) {  // 得到 2根
        double lnfpg = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
        double lnfpl = Zl - 1 - log(Zl - B_) - A_ / (2 * sqrt(2) * B_) * log((Zl + B_ + sqrt(2) * B_) / (Zl + B_ - sqrt(2) * B_));
        Z_ = (lnfpg - lnfpl < 0) ? Zg : Zl;  // 保留逸度较小的结果 可能报错
    } else {
        Z_ = Zg;
    }
}

void PRimp::Tflash(double Ts) {
    is_single_phase_ = false;  // Tflash 得到饱和线
    T_ = Ts;
    double a = ac_ * pow(1 + kappa_ * (1 - sqrt(T_ / Tc_)), 2);
    double b = bc_;
    double lnfpg, lnfpl, lnfp, lnfp_pmin, lnfp_pmax;  // 记录参数
    double ps, ps_min, ps_max;  // 二分法计算
    double Z, Zg, Zl;  // 用于储存一元三次方程求解结果
    // 1.寻找下边界
    ps_min = 1;  // 1Pa 开始
    A_ = a * ps_min / pow(R_, 2) / pow(T_, 2);
    B_ = b * ps_min / R_ / T_;
    if (ShengjinFanSolver(1, -(1 - B_), A_ - 3 * pow(B_, 2) - 2 * B_, -(A_ * B_ - pow(B_, 2) - pow(B_, 3)), &Zg, &Zl) == true) {  // 得到 2根
        lnfpg = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
        lnfpl = Zl - 1 - log(Zl - B_) - A_ / (2 * sqrt(2) * B_) * log((Zl + B_ + sqrt(2) * B_) / (Zl + B_ - sqrt(2) * B_));
        lnfp = lnfpg - lnfpl;
    } else {
        lnfp = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
        if (Zg < Zc_) lnfp = -lnfp;  // 液相则异号
    }
    lnfp_pmin = lnfp;
    // 2.寻找上边界
    ps_max = pc_ - 1;  // 低于临界温度 1Pa
    A_ = a * ps_max / pow(R_, 2) / pow(T_, 2);
    B_ = b * ps_max / R_ / T_;
    if (ShengjinFanSolver(1, -(1 - B_), A_ - 3 * pow(B_, 2) - 2 * B_, -(A_ * B_ - pow(B_, 2) - pow(B_, 3)), &Zg, &Zl) == true) {  // 得到 2根
        lnfpg = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
        lnfpl = Zl - 1 - log(Zl - B_) - A_ / (2 * sqrt(2) * B_) * log((Zl + B_ + sqrt(2) * B_) / (Zl + B_ - sqrt(2) * B_));
        lnfp = lnfpg - lnfpl;
    } else {
        lnfp = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
        if (Zg < Zc_) lnfp = -lnfp;  // 液相则异号
    }
    lnfp_pmax = lnfp;
    // 3.判断是否有解
    if (lnfp_pmin * lnfp_pmax > 0 && T_ < 0.95 * Tc_) {
        T_ = 0.6 * Tc_;
        P_ = 0.001;
        Zg_ = 1.0;
        Zl_ = 1E-9;
        return;  // 计算失败 返回近零压力点
    }
    if (lnfp_pmin * lnfp_pmax > 0 && T_ > 0.95 * Tc_) {
        T_ = Tc_;
        P_ = pc_;
        Zg_ = Zc_;
        Zl_ = Zc_;
        return;  // 计算失败 返回临界点
    }
    // 4.二分法求解：
    // 随着压力的增加 lnfp 单调增加 就是这么离谱 不要问我为什么 因为我也不理解
    int counter = 0;
    while (1) {
        if (counter == 1E3) {
            T_ = Tc_;
            P_ = pc_;
            Zg_ = Zc_;
            Zl_ = Zc_;
            return;  // 循环次数过多 返回临界点
        }
        counter++;
        ps = (ps_min + ps_max) / 2;
        A_ = a * ps / pow(R_, 2) / pow(T_, 2);
        B_ = b * ps / R_ / T_;
        if (ShengjinFanSolver(1, -(1 - B_), A_ - 3 * pow(B_, 2) - 2 * B_, -(A_ * B_ - pow(B_, 2) - pow(B_, 3)), &Zg, &Zl) == true) {  // 得到 2根
            lnfpg = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
            lnfpl = Zl - 1 - log(Zl - B_) - A_ / (2 * sqrt(2) * B_) * log((Zl + B_ + sqrt(2) * B_) / (Zl + B_ - sqrt(2) * B_));
            lnfp = lnfpg - lnfpl;
        } else {
            lnfp = Zg - 1 - log(Zg - B_) - A_ / (2 * sqrt(2) * B_) * log((Zg + B_ + sqrt(2) * B_) / (Zg + B_ - sqrt(2) * B_));
            if (Zg < Zc_) lnfp = -lnfp;  // 液相则异号
        }
        if (abs(lnfp) < 1E-6) {
            P_ = ps;
            Zg_ = Zg;
            Zl_ = Zl;
            return;  // 计算收敛 返回结果
        } else if (lnfp * lnfp_pmin < 0) {
            lnfp_pmax = lnfp;
            ps_max = ps;
        } else if (lnfp * lnfp_pmax < 0) {
            lnfp_pmin = lnfp;
            ps_min = ps;
        }
    }
}

double PRimp::Temperature() {
    return T_;
}

double PRimp::Pressure() {
    return is_single_phase_ ? P_ : 0;
}

double PRimp::Density() {
    return is_single_phase_ ? P_ / (Z_ * R_ * T_) : 0;
}

double PRimp::Ps() {
    return is_single_phase_ ? 0 : P_;
}

double PRimp::Rhog() {
    return is_single_phase_ ? 0 : P_ / (Zg_ * R_ * T_);
}

double PRimp::Rhol() {
    return is_single_phase_ ? 0 : P_ / (Zl_ * R_ * T_);
}

bool PRimp::ShengjinFanSolver(double a, double b, double c, double d,
                              double* Zg, double* Zl) {
// 返回值 bool has_double_root
    double A = pow(b, 2) - 3 * a * c;
    double B = b * c - 9 * a * d;
    double C = pow(c, 2) - 3 * b * d;
    double Delta = pow(B, 2) - 4 * A * C;
    double Z = 0;
    if (abs(A) < 1E-16 && abs(B) < 1E-16) {
        // A=B=0 存在 三重实根 临界点？
        Z = -b / (3 * a);
    } else if (Delta > 0) {
        // Delta>0 方程有一个实根和一对共轭复根 仅保留实根！!
        double Y1 = A * b + 1.5 * a * (-B + sqrt(Delta));
        double Y2 = A * b + 1.5 * a * (-B - sqrt(Delta));
        Z = (-b - (cbrt(Y1) + cbrt(Y2))) / (3 * a);
    } else if (abs(Delta) < 1E-16) {
        // Delta=0 方程有三个实根 其中有一个二重根 舍弃二重根！！
        Z = B / A - b / a;;
    } else if (Delta < 0) {
        // Delta<0 方程有三个不相等的实根 舍弃中间根！！
        double theta3 = acos((2 * A * b - 3 * a * B) / (2 * A * sqrt(A))) / 3;
        double x1 = (-b - 2 * sqrt(A) * cos(theta3)) / (3 * a);
        double x2 = (-b + sqrt(A) * (cos(theta3) + sqrt(3) * sin(theta3))) / (3 * a);
        double x3 = (-b + sqrt(A) * (cos(theta3) - sqrt(3) * sin(theta3))) / (3 * a);
        double Zg0 = (x1 > x2) ? (x1 > x3 ? x1 : x3) : (x2 > x3 ? x2 : x3);
        double Zl0 = (x1 < x2) ? (x1 < x3 ? x1 : x3) : (x2 < x3 ? x2 : x3);
        if (Zl0 < 0) {
            Z = Zg0;
        } else {
            *Zg = Zg0;
            *Zl = Zl0;
            return true;
        }
    }
    *Zg = Z;
    *Zl = Z;
    return false;
}

}  // namespace thermolib

// 以下内容用于第三方封装：c++
extern "C" __declspec(dllexport) PR * __stdcall GetPR() {
    return new thermolib::PRimp;
}

// 以下内容用于第三方封装：python
extern "C" {
    thermolib::PRimp pr;
    __declspec(dllexport) void SetFluid(double Tc, double pc, double R, double omega) {
        return pr.SetFluid(Tc, pc, R, omega);
    }  // 加载方程参数
    __declspec(dllexport) void TPflash(double T, double p) { return pr.TPflash(T, p); }
    __declspec(dllexport) void Tflash(double Ts) { return pr.Tflash(Ts); }
    __declspec(dllexport) double Temperature() { return pr.Temperature(); }
    __declspec(dllexport) double Pressure() { return pr.Pressure(); }
    __declspec(dllexport) double Density() { return pr.Density(); }
    __declspec(dllexport) double Ps() { return pr.Ps(); }
    __declspec(dllexport) double Rhog() { return pr.Rhog(); }
    __declspec(dllexport) double Rhol() { return pr.Rhol(); }
}

