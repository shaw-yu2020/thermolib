
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

#ifndef SRC_RK_RK_IMP_H_
#define SRC_RK_RK_IMP_H_

#include<cmath>
#include"src/rk/rk.h"

namespace thermolib {

using std::abs;

class RKimp:public RK {
public:
    void SetFluid(double Tc, double pc, double R) override;
    void TPflash(double T, double p) override;
    void Tflash(double Ts) override;
    double Temperature() override;
    double Pressure() override;
    double Density() override;
    double Ps() override;
    double Rhog() override;
    double Rhol() override;
private:
    bool is_single_phase_ = true;
    double A_, B_, Z_, T_, P_, Zg_, Zl_;  // 记录参数
    double Tc_, pc_, R_, a_, b_, Zc_;
    bool ShengjinFanSolver(double a, double b, double c, double d, double* Zg, double* Zl);  // 返回值 bool has_double_root
};

}  // namespace thermolib

#endif  // SRC_RK_RK_IMP_H_

