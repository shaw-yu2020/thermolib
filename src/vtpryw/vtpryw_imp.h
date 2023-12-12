
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

#ifndef SRC_VTPRYW_VTPRYW_IMP_H_
#define SRC_VTPRYW_VTPRYW_IMP_H_

#include<cmath>
#include"src/vtpryw/vtpryw.h"

namespace thermolib {

using std::abs;

class VTPRYWimp:public VTPRYW {
public:
    void SetFluid(double Tc, double pc, double R, double omega,
                  double Z, double beta, double gamma, double eta, double xi) override;
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
    double A_, B_, C_, Z_, T_, P_, Zg_, Zl_;
    double Tc_, pc_, R_, kappa_, ac_, bc_, cc_, Zc_, beta_, gamma_, eta_, xi_;
    bool ShengjinFanSolver(double a, double b, double c, double d, double* Zg, double* Zl);  // 返回值 bool has_double_root
};

}  // namespace thermolib

#endif  // SRC_VTPRYW_VTPRYW_IMP_H_

