
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

#ifndef SRC_PR_PR_IMP_H_
#define SRC_PR_PR_IMP_H_

#include<cmath>
#include"src/pr/pr.h"

namespace thermolib {

using std::abs;

class PRimp:public PR {
public:
    void SetFluid(double Tc, double pc, double R, double omega) override;
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
    double Tc_ = 0, pc_ = 0, R_, kappa_ = 0, ac_ = 0, bc_ = 0, Zc_ = 0;  // omega 仅仅用于计算 kappa
    bool ShengjinFanSolver(double a, double b, double c, double d, double* Zg, double* Zl);  // 返回值 bool has_double_root
};

}  // namespace thermolib

#endif  // SRC_PR_PR_IMP_H_

