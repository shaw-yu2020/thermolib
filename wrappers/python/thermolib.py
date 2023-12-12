'''
Copyright[2022]
Author: ygl
Email: ygl_2020@stu.xjtu.edu.cn
'''

import ctypes

def get_rk(path):
    rk_dll = ctypes.CDLL(path, winmode=0)
    rk_dll.SetFluid.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double)
    rk_dll.TPflash.argtypes = (ctypes.c_double, ctypes.c_double)
    rk_dll.Tflash.argtypes = (ctypes.c_double, )
    rk_dll.Temperature.restype = ctypes.c_double
    rk_dll.Pressure.restype = ctypes.c_double
    rk_dll.Density.restype = ctypes.c_double
    rk_dll.Ps.restype = ctypes.c_double
    rk_dll.Rhog.restype = ctypes.c_double
    rk_dll.Rhol.restype = ctypes.c_double
    return rk_dll

def get_srk(path):
    srk_dll = ctypes.CDLL(path, winmode=0)
    srk_dll.SetFluid.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    srk_dll.TPflash.argtypes = (ctypes.c_double, ctypes.c_double)
    srk_dll.Tflash.argtypes = (ctypes.c_double, )
    srk_dll.Temperature.restype = ctypes.c_double
    srk_dll.Pressure.restype = ctypes.c_double
    srk_dll.Density.restype = ctypes.c_double
    srk_dll.Ps.restype = ctypes.c_double
    srk_dll.Rhog.restype = ctypes.c_double
    srk_dll.Rhol.restype = ctypes.c_double
    return srk_dll

def get_pr(path):
    pr_dll = ctypes.CDLL(path, winmode=0)
    pr_dll.SetFluid.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    pr_dll.TPflash.argtypes = (ctypes.c_double, ctypes.c_double)
    pr_dll.Tflash.argtypes = (ctypes.c_double, )
    pr_dll.Temperature.restype = ctypes.c_double
    pr_dll.Pressure.restype = ctypes.c_double
    pr_dll.Density.restype = ctypes.c_double
    pr_dll.Ps.restype = ctypes.c_double
    pr_dll.Rhog.restype = ctypes.c_double
    pr_dll.Rhol.restype = ctypes.c_double
    return pr_dll

def main():
    return 0

if __name__ == '__main__':
    main()

