
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

#include"src/srk/srk_imp.h"

namespace thermolib {

void SRKimp::SetFluid(double Tc, double pc, double R, double omega) {
    Tc_ = Tc;
    pc_ = pc;
    R_ = R;
    m_ = 0.480 + 1.574 * omega - 0.176 * pow(omega, 2);
    ac_ = 0.42747 * pow(R, 2) * pow(Tc, 2) / pc;
    bc_ = 0.08664 * R * Tc / pc;
    Zc_ = 1.0 / 3;
}

void SRKimp::TPflash(double T, double p) {
    is_single_phase_ = true;  // TPflash 得到单相区
    T_ = T;
    P_ = p;
    double a = ac_ * pow(1 + m_ * (1 - sqrt(T / Tc_)), 2);
    double b = bc_;
    A_ = a * P_ / pow(R_, 2) / pow(T_, 2);
    B_ = b * P_ / R_ / T_;
    double Zg = 0, Zl = 0;  // 用于记录一元三次方程求解结果
    if (ShengjinFanSolver(1, -1, A_ - B_ - pow(B_, 2), -A_ * B_, &Zg, &Zl) == true) {  // 得到 2根
        double lnfpg = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
        double lnfpl = Zl - 1 - log(Zl - B_) - A_ / B_ * log(Zl / (Zl + B_));
        Z_ = (lnfpg - lnfpl < 0) ? Zg : Zl;  // 保留逸度较小的结果 可能报错
    } else {
        Z_ = Zg;
    }
}

void SRKimp::Tflash(double Ts) {
    is_single_phase_ = false;  // Tflash 得到饱和线
    T_ = Ts;
    double a = ac_ * pow(1 + m_ * (1 - sqrt(Ts / Tc_)), 2);
    double b = bc_;
    double lnfpg, lnfpl, lnfp, lnfp_pmin, lnfp_pmax;  // 记录参数
    double ps, ps_min, ps_max;  // 二分法计算
    double Z, Zg, Zl;  // 用于储存一元三次方程求解结果
    // 1.寻找下边界
    ps_min = 1;  // 1Pa 开始
    A_ = a * ps_min / pow(R_, 2) / pow(T_, 2.5);
    B_ = b * ps_min / R_ / T_;
    if (ShengjinFanSolver(1, -1, A_ - B_ - pow(B_, 2), -A_ * B_, &Zg, &Zl) == true) {  // 得到 2根
        lnfpg = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
        lnfpl = Zl - 1 - log(Zl - B_) - A_ / B_ * log(Zl / (Zl + B_));
        lnfp = lnfpg - lnfpl;
    } else {
        lnfp = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
        if (Zg < Zc_) lnfp = -lnfp;  // 液相则异号
    }
    lnfp_pmin = lnfp;
    // 2.寻找上边界
    ps_max = pc_ - 1;  // 低于临界温度 1Pa
    A_ = a * ps_max / pow(R_, 2) / pow(T_, 2.5);
    B_ = b * ps_max / R_ / T_;
    if (ShengjinFanSolver(1, -1, A_ - B_ - pow(B_, 2), -A_ * B_, &Zg, &Zl) == true) {  // 得到 2根
        lnfpg = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
        lnfpl = Zl - 1 - log(Zl - B_) - A_ / B_ * log(Zl / (Zl + B_));
        lnfp = lnfpg - lnfpl;
    } else {
        lnfp = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
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
        A_ = a * ps / pow(R_, 2) / pow(T_, 2.5);
        B_ = b * ps / R_ / T_;
        if (ShengjinFanSolver(1, -1, A_ - B_ - pow(B_, 2), -A_ * B_, &Zg, &Zl) == true) {  // 得到 2根
            lnfpg = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
            lnfpl = Zl - 1 - log(Zl - B_) - A_ / B_ * log(Zl / (Zl + B_));
            lnfp = lnfpg - lnfpl;
        } else {
            lnfp = Zg - 1 - log(Zg - B_) - A_ / B_ * log(Zg / (Zg + B_));
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

double SRKimp::Temperature() {
    return T_;
}

double SRKimp::Pressure() {
    return is_single_phase_ ? P_ : 0;
}

double SRKimp::Density() {
    return is_single_phase_ ? P_ / (Z_ * R_ * T_) : 0;
}

double SRKimp::Ps() {
    return is_single_phase_ ? 0 : P_;
}

double SRKimp::Rhog() {
    return is_single_phase_ ? 0 : P_ / (Zg_ * R_ * T_);
}

double SRKimp::Rhol() {
    return is_single_phase_ ? 0 : P_ / (Zl_ * R_ * T_);
}

bool SRKimp::ShengjinFanSolver(double a, double b, double c, double d,
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
extern "C" __declspec(dllexport) SRK * __stdcall GetSRK() {
    return new thermolib::SRKimp;
}

// 以下内容用于第三方封装：python
extern "C" {
    thermolib::SRKimp srk;
    __declspec(dllexport) void SetFluid(double Tc, double pc, double R, double omega) {
        return srk.SetFluid(Tc, pc, R, omega);
    }  // 加载方程参数
    __declspec(dllexport) void TPflash(double T, double p) { return srk.TPflash(T, p); }
    __declspec(dllexport) void Tflash(double Ts) { return srk.Tflash(Ts); }
    __declspec(dllexport) double Temperature() { return srk.Temperature(); }
    __declspec(dllexport) double Pressure() { return srk.Pressure(); }
    __declspec(dllexport) double Density() { return srk.Density(); }
    __declspec(dllexport) double Ps() { return srk.Ps(); }
    __declspec(dllexport) double Rhog() { return srk.Rhog(); }
    __declspec(dllexport) double Rhol() { return srk.Rhol(); }
}

