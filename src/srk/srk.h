
// Copyright[2022]
// Author: ygl
// Email: ygl_2020@stu.xjtu.edu.cn

// 避免多次访问头文件
#ifndef SRC_SRK_SRK_H_
#define SRC_SRK_SRK_H_

class SRK {
public:
    virtual void SetFluid(double Tc, double pc, double R, double omega) = 0;
    virtual void TPflash(double T, double p) = 0;
    virtual void Tflash(double Ts) = 0;
    virtual double Temperature() = 0;
    virtual double Pressure() = 0;
    virtual double Density() = 0;
    virtual double Ps() = 0;
    virtual double Rhog() = 0;
    virtual double Rhol() = 0;
};

extern "C" __declspec(dllexport) SRK * __stdcall GetSRK();

#endif  // SRC_SRK_SRK_H_

