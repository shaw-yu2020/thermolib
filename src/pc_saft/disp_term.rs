use std::f64::consts::{PI, TAU};
pub struct DispTerm {
    m: f64,
    c: CTerm,   // CTerm
    i1: I1Term, // I1Term
    i2: I2Term, // I2Term
    // cached variables
    c1t0d1: (f64, f64),
    c1t0d2: (f64, f64),
    c1t0d3: (f64, f64),
    c1t0d4: (f64, f64),
    c1t1d0: (f64, f64),
    c1t1d1: (f64, f64),
    c1t1d2: (f64, f64),
    c1t1d3: (f64, f64),
}
impl DispTerm {
    pub fn new(m: f64) -> Self {
        Self {
            m,
            c: CTerm::new(m),   // CTerm
            i1: I1Term::new(m), // I1Term
            i2: I2Term::new(m), // I1Term
            // cached variables
            c1t0d1: (0.0, 0.0),
            c1t0d2: (0.0, 0.0),
            c1t0d3: (0.0, 0.0),
            c1t0d4: (0.0, 0.0),
            c1t1d0: (0.0, 0.0),
            c1t1d1: (0.0, 0.0),
            c1t1d2: (0.0, 0.0),
            c1t1d3: (0.0, 0.0),
        }
    }
    pub fn t0d0(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64) -> f64 {
        -TAU * m2s3e1 * self.i1.t0d0(eta)
            - PI * m2s3e2 * self.m * self.c1t0d0(eta) * self.i2.t0d0(eta)
    }
    pub fn t0d1(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64) -> f64 {
        -TAU * m2s3e1 * self.i1.t0d1(eta)
            - (PI * m2s3e2 * self.m)
                * (self.c1t0d1(eta) * self.i2.t0d0(eta) + self.c1t0d0(eta) * self.i2.t0d1(eta))
    }
    pub fn t0d2(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64) -> f64 {
        -TAU * m2s3e1 * self.i1.t0d2(eta)
            - (PI * m2s3e2 * self.m)
                * (self.c1t0d2(eta) * self.i2.t0d0(eta)
                    + 2.0 * self.c1t0d1(eta) * self.i2.t0d1(eta)
                    + self.c1t0d0(eta) * self.i2.t0d2(eta))
    }
    pub fn t0d3(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64) -> f64 {
        -TAU * m2s3e1 * self.i1.t0d3(eta)
            - (PI * m2s3e2 * self.m)
                * (self.c1t0d3(eta) * self.i2.t0d0(eta)
                    + 3.0 * self.c1t0d2(eta) * self.i2.t0d1(eta)
                    + 3.0 * self.c1t0d1(eta) * self.i2.t0d2(eta)
                    + self.c1t0d0(eta) * self.i2.t0d3(eta))
    }
    pub fn t0d4(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64) -> f64 {
        -TAU * m2s3e1 * self.i1.t0d4(eta)
            - (PI * m2s3e2 * self.m)
                * (self.c1t0d4(eta) * self.i2.t0d0(eta)
                    + 4.0 * self.c1t0d3(eta) * self.i2.t0d1(eta)
                    + 6.0 * self.c1t0d2(eta) * self.i2.t0d2(eta)
                    + 4.0 * self.c1t0d1(eta) * self.i2.t0d3(eta)
                    + self.c1t0d0(eta) * self.i2.t0d4(eta))
    }
    pub fn t1d0(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64, eta1: f64) -> f64 {
        -TAU * m2s3e1 * (eta1 * self.i1.t1d0(eta) - self.i1.t0d0(eta))
            - (PI * m2s3e2 * self.m)
                * (eta1
                    * (self.c1t1d0(eta) * self.i2.t0d0(eta) + self.c1t0d0(eta) * self.i2.t1d0(eta))
                    - 2.0 * self.c1t0d0(eta) * self.i2.t0d0(eta))
    }
    pub fn t1d1(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64, eta1: f64) -> f64 {
        -TAU * m2s3e1 * (eta1 * self.i1.t1d1(eta) - self.i1.t0d1(eta))
            - PI * m2s3e2
                * self.m
                * (eta1
                    * (self.c1t1d1(eta) * self.i2.t0d0(eta)
                        + self.c1t1d0(eta) * self.i2.t0d1(eta)
                        + self.c1t0d1(eta) * self.i2.t1d0(eta)
                        + self.c1t0d0(eta) * self.i2.t1d1(eta))
                    - 2.0
                        * (self.c1t0d1(eta) * self.i2.t0d0(eta)
                            + self.c1t0d0(eta) * self.i2.t0d1(eta)))
    }
    pub fn t1d2(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64, eta1: f64) -> f64 {
        -TAU * m2s3e1 * (eta1 * self.i1.t1d2(eta) - self.i1.t0d2(eta))
            - (PI * m2s3e2 * self.m)
                * (eta1
                    * (self.c1t1d2(eta) * self.i2.t0d0(eta)
                        + 2.0 * self.c1t1d1(eta) * self.i2.t0d1(eta)
                        + self.c1t1d0(eta) * self.i2.t0d2(eta)
                        + self.c1t0d2(eta) * self.i2.t1d0(eta)
                        + 2.0 * self.c1t0d1(eta) * self.i2.t1d1(eta)
                        + self.c1t0d0(eta) * self.i2.t1d2(eta))
                    - 2.0
                        * (self.c1t0d2(eta) * self.i2.t0d0(eta)
                            + 2.0 * self.c1t0d1(eta) * self.i2.t0d1(eta)
                            + self.c1t0d0(eta) * self.i2.t0d2(eta)))
    }
    pub fn t1d3(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64, eta1: f64) -> f64 {
        -TAU * m2s3e1 * (eta1 * self.i1.t1d3(eta) - self.i1.t0d3(eta))
            - (PI * m2s3e2 * self.m)
                * (eta1
                    * (self.c1t1d3(eta) * self.i2.t0d0(eta)
                        + 3.0 * self.c1t1d2(eta) * self.i2.t0d1(eta)
                        + 3.0 * self.c1t1d1(eta) * self.i2.t0d2(eta)
                        + self.c1t1d0(eta) * self.i2.t0d3(eta)
                        + self.c1t0d3(eta) * self.i2.t1d0(eta)
                        + 3.0 * self.c1t0d2(eta) * self.i2.t1d1(eta)
                        + 3.0 * self.c1t0d1(eta) * self.i2.t1d2(eta)
                        + self.c1t0d0(eta) * self.i2.t1d3(eta))
                    - 2.0
                        * (self.c1t0d3(eta) * self.i2.t0d0(eta)
                            + 3.0 * self.c1t0d2(eta) * self.i2.t0d1(eta)
                            + 3.0 * self.c1t0d1(eta) * self.i2.t0d2(eta)
                            + self.c1t0d0(eta) * self.i2.t0d3(eta)))
    }
    pub fn t2d0(&mut self, m2s3e1: f64, m2s3e2: f64, eta: f64, eta1: f64, eta2: f64) -> f64 {
        -TAU * m2s3e1
            * (self.i1.t2d0(eta, eta1, eta2) - 2.0 * eta1 * self.i1.t1d0(eta)
                + 2.0 * self.i1.t0d0(eta))
            - (PI * m2s3e2 * self.m)
                * (self.c1t2d0(eta, eta1, eta2) * self.i2.t0d0(eta)
                    + 2.0 * eta1.powi(2) * self.c1t1d0(eta) * self.i2.t1d0(eta)
                    + self.c1t0d0(eta) * self.i2.t2d0(eta, eta1, eta2)
                    - (4.0 * eta1)
                        * (self.c1t1d0(eta) * self.i2.t0d0(eta)
                            + self.c1t0d0(eta) * self.i2.t1d0(eta))
                    + 6.0 * self.c1t0d0(eta) * self.i2.t0d0(eta))
    }
}
impl DispTerm {
    #[inline]
    fn c1t0d0(&mut self, eta: f64) -> f64 {
        self.c.eta0(eta).recip()
    }
    fn c1t0d1(&mut self, eta: f64) -> f64 {
        if self.c1t0d1.0 != eta {
            self.c1t0d1 = (eta, -eta / self.c.eta0(eta).powi(2) * self.c.eta1(eta))
        }
        self.c1t0d1.1
    }
    fn c1t0d2(&mut self, eta: f64) -> f64 {
        if self.c1t0d2.0 != eta {
            self.c1t0d2 = (
                eta,
                -eta.powi(2) / self.c.eta0(eta).powi(3)
                    * (self.c.eta2(eta) * self.c.eta0(eta) - 2.0 * self.c.eta1(eta).powi(2)),
            )
        }
        self.c1t0d2.1
    }
    fn c1t0d3(&mut self, eta: f64) -> f64 {
        if self.c1t0d3.0 != eta {
            self.c1t0d3 = (
                eta,
                -eta.powi(3) / self.c.eta0(eta).powi(4)
                    * (self.c.eta3(eta) * self.c.eta0(eta).powi(2)
                        - 6.0 * self.c.eta2(eta) * self.c.eta1(eta) * self.c.eta0(eta)
                        + 6.0 * self.c.eta1(eta).powi(3)),
            )
        }
        self.c1t0d3.1
    }
    fn c1t0d4(&mut self, eta: f64) -> f64 {
        if self.c1t0d4.0 != eta {
            self.c1t0d4 = (
                eta,
                -eta.powi(4) / self.c.eta0(eta).powi(5)
                    * (self.c.eta4(eta) * self.c.eta0(eta).powi(3)
                        - 8.0 * self.c.eta3(eta) * self.c.eta1(eta) * self.c.eta0(eta).powi(2)
                        - 6.0 * self.c.eta2(eta).powi(2) * self.c.eta0(eta).powi(2)
                        + 36.0 * self.c.eta2(eta) * self.c.eta1(eta).powi(2) * self.c.eta0(eta)
                        - 24.0 * self.c.eta1(eta).powi(4)),
            )
        }
        self.c1t0d4.1
    }
    fn c1t1d0(&mut self, eta: f64) -> f64 {
        if self.c1t1d0.0 != eta {
            self.c1t1d0 = (eta, -self.c.eta0(eta).powi(2).recip() * self.c.eta1(eta))
        }
        self.c1t1d0.1
    }
    fn c1t1d1(&mut self, eta: f64) -> f64 {
        if self.c1t1d1.0 != eta {
            self.c1t1d1 = (
                eta,
                -self.c.eta0(eta).powi(3).recip()
                    * (eta
                        * (self.c.eta2(eta) * self.c.eta0(eta) - 2.0 * self.c.eta1(eta).powi(2))
                        + self.c.eta1(eta) * self.c.eta0(eta)),
            )
        }
        self.c1t1d1.1
    }
    fn c1t1d2(&mut self, eta: f64) -> f64 {
        if self.c1t1d2.0 != eta {
            self.c1t1d2 = (
                eta,
                -eta / self.c.eta0(eta).powi(4)
                    * (eta
                        * (self.c.eta3(eta) * self.c.eta0(eta).powi(2)
                            - 6.0 * self.c.eta2(eta) * self.c.eta1(eta) * self.c.eta0(eta)
                            + 6.0 * self.c.eta1(eta).powi(3))
                        + 2.0 * self.c.eta2(eta) * self.c.eta0(eta).powi(2)
                        - 4.0 * self.c.eta1(eta).powi(2) * self.c.eta0(eta)),
            )
        }
        self.c1t1d2.1
    }
    fn c1t1d3(&mut self, eta: f64) -> f64 {
        if self.c1t1d3.0 != eta {
            self.c1t1d3 = (
                eta,
                -eta.powi(2) / self.c.eta0(eta).powi(5)
                    * (eta
                        * (self.c.eta0(eta).powi(3) * self.c.eta4(eta)
                            - self.c.eta0(eta).powi(2)
                                * (8.0 * self.c.eta3(eta) * self.c.eta1(eta)
                                    + 6.0 * self.c.eta2(eta).powi(2))
                            + self.c.eta0(eta)
                                * (36.0 * self.c.eta2(eta) * self.c.eta1(eta).powi(2))
                            - 24.0 * self.c.eta1(eta).powi(4))
                        + 3.0 * self.c.eta3(eta) * self.c.eta0(eta).powi(3)
                        - 18.0 * self.c.eta2(eta) * self.c.eta1(eta) * self.c.eta0(eta).powi(2)
                        + 18.0 * self.c.eta1(eta) * self.c.eta1(eta).powi(2) * self.c.eta0(eta)),
            )
        }
        self.c1t1d3.1
    }
    fn c1t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        2.0 * (eta1 * self.c.eta1(eta)).powi(2) / self.c.eta0(eta).powi(3)
            - (eta2 * self.c.eta1(eta) + eta1.powi(2) * self.c.eta2(eta)) / self.c.eta0(eta).powi(2)
    }
}
/// CTerm
struct CTerm {
    m: f64,
    // cached variables
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    eta3: (f64, f64),
    eta4: (f64, f64),
}
impl CTerm {
    fn new(m: f64) -> Self {
        Self {
            m,
            // cached variables
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            eta3: (0.0, 0.0),
            eta4: (0.0, 0.0),
        }
    }
    fn eta0(&mut self, eta: f64) -> f64 {
        if eta != self.eta0.0 {
            self.eta0 = (
                eta,
                1.0 + 2.0 * self.m * (4.0 * eta - eta.powi(2)) / (1.0 - eta).powi(4)
                    + (1.0 - self.m)
                        * ((20.0 * eta + 12.0 * eta.powi(3))
                            - (27.0 * eta.powi(2) + 2.0 * eta.powi(4)))
                        / ((1.0 - eta) * (2.0 - eta)).powi(2),
            )
        }
        self.eta0.1
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        if eta != self.eta1.0 {
            self.eta1 = (
                eta,
                4.0 * self.m * (-eta.powi(2) + 5.0 * eta + 2.0) / (1.0 - eta).powi(5)
                    + 2.0 * (1.0 - self.m) * (eta.powi(3) + 6.0 * eta.powi(2) - 24.0 * eta + 20.0)
                        / ((1.0 - eta) * (2.0 - eta)).powi(3),
            )
        }
        self.eta1.1
    }
    fn eta2(&mut self, eta: f64) -> f64 {
        if eta != self.eta2.0 {
            self.eta2 = (
                eta,
                12.0 * self.m * (-eta.powi(2) + 6.0 * eta + 5.0) / (1.0 - eta).powi(6)
                    + (6.0 * (1.0 - self.m))
                        * ((48.0 * eta.powi(2) + 44.0)
                            - (eta.powi(4) + 8.0 * eta.powi(3) + 80.0 * eta))
                        / ((1.0 - eta) * (2.0 - eta)).powi(4),
            )
        }
        self.eta2.1
    }
    fn eta3(&mut self, eta: f64) -> f64 {
        if eta != self.eta3.0 {
            self.eta3 = (
                eta,
                48.0 * self.m * (-eta.powi(2) + 7.0 * eta + 9.0) / (1.0 - eta).powi(7)
                    + (24.0 * (1.0 - self.m))
                        * ((eta.powi(5) + 10.0 * eta.powi(4) + 200.0 * eta.powi(2) + 92.0)
                            - (80.0 * eta.powi(3) + 220.0 * eta))
                        / ((1.0 - eta) * (2.0 - eta)).powi(5),
            )
        }
        self.eta3.1
    }
    fn eta4(&mut self, eta: f64) -> f64 {
        if eta != self.eta4.0 {
            self.eta4 = (
                eta,
                240.0 * self.m * (-eta.powi(2) + 8.0 * eta + 14.0) / (1.0 - eta).powi(8)
                    + (120.0 * (1.0 - self.m))
                        * ((120.0 * eta.powi(4) + 660.0 * eta.powi(2) + 188.0)
                            - (eta.powi(6)
                                + 12.0 * eta.powi(5)
                                + 400.0 * eta.powi(3)
                                + 552.0 * eta))
                        / ((1.0 - eta) * (2.0 - eta)).powi(6),
            )
        }
        self.eta4.1
    }
}
/// I1Term
struct I1Term {
    a0: f64,
    a1: f64,
    a2: f64,
    a3: f64,
    a4: f64,
    a5: f64,
    a6: f64,
    // cached variables 0
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    // cached variables 1
    t0d1: (f64, f64),
    t0d2: (f64, f64),
    t0d3: (f64, f64),
    t0d4: (f64, f64),
    t1d1: (f64, f64),
    t1d2: (f64, f64),
    t1d3: (f64, f64),
}
// Universal Model Constants
const A00: f64 = 0.9105631445;
const A01: f64 = 0.6361281449;
const A02: f64 = 2.6861347891;
const A03: f64 = -26.547362491;
const A04: f64 = 97.759208784;
const A05: f64 = -159.59154087;
const A06: f64 = 91.297774084;
const A10: f64 = -0.3084016918;
const A11: f64 = 0.1860531159;
const A12: f64 = -2.5030047259;
const A13: f64 = 21.419793629;
const A14: f64 = -65.255885330;
const A15: f64 = 83.318680481;
const A16: f64 = -33.746922930;
const A20: f64 = -0.0906148351;
const A21: f64 = 0.4527842806;
const A22: f64 = 0.5962700728;
const A23: f64 = -1.7241829131;
const A24: f64 = -4.1302112531;
const A25: f64 = 13.776631870;
const A26: f64 = -8.6728470368;
impl I1Term {
    fn new(m: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            a0: A00 + m1 * A10 + m12 * A20,
            a1: A01 + m1 * A11 + m12 * A21,
            a2: A02 + m1 * A12 + m12 * A22,
            a3: A03 + m1 * A13 + m12 * A23,
            a4: A04 + m1 * A14 + m12 * A24,
            a5: A05 + m1 * A15 + m12 * A25,
            a6: A06 + m1 * A16 + m12 * A26,
            // cached variables 0
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            // cached variables 1
            t0d1: (0.0, 0.0),
            t0d2: (0.0, 0.0),
            t0d3: (0.0, 0.0),
            t0d4: (0.0, 0.0),
            t1d1: (0.0, 0.0),
            t1d2: (0.0, 0.0),
            t1d3: (0.0, 0.0),
        }
    }
    fn eta0(&mut self, eta: f64) -> f64 {
        if eta != self.eta0.0 {
            self.eta0 = (
                eta,
                self.a0
                    + self.a1 * eta
                    + self.a2 * eta.powi(2)
                    + self.a3 * eta.powi(3)
                    + self.a4 * eta.powi(4)
                    + self.a5 * eta.powi(5)
                    + self.a6 * eta.powi(6),
            )
        }
        self.eta0.1
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        if eta != self.eta1.0 {
            self.eta1 = (
                eta,
                self.a1
                    + self.a2 * 2.0 * eta
                    + self.a3 * 3.0 * eta.powi(2)
                    + self.a4 * 4.0 * eta.powi(3)
                    + self.a5 * 5.0 * eta.powi(4)
                    + self.a6 * 6.0 * eta.powi(5),
            )
        }
        self.eta1.1
    }
    fn eta2(&mut self, eta: f64) -> f64 {
        if eta != self.eta2.0 {
            self.eta2 = (
                eta,
                self.a2 * 2.0
                    + self.a3 * 6.0 * eta
                    + self.a4 * 12.0 * eta.powi(2)
                    + self.a5 * 20.0 * eta.powi(3)
                    + self.a6 * 30.0 * eta.powi(4),
            )
        }
        self.eta2.1
    }
    /// equal to = [rho*i1]_t0d0 / rho
    /// equal to = i1t0d0
    /// equal to = eta0
    #[inline]
    fn t0d0(&mut self, eta: f64) -> f64 {
        self.eta0(eta)
    }
    /// equal to = { [rho*i1]_t1d0 / rho } / eta1
    /// equal to = { i1t1d0 } / eta1
    /// equal to = eta1
    #[inline]
    fn t1d0(&mut self, eta: f64) -> f64 {
        self.eta1(eta)
    }
    /// equal to = [rho*i1]_t2d0 / rho
    /// equal to = i1t2d0
    fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        eta2 * self.eta1(eta) + eta1.powi(2) * self.eta2(eta)
    }
    /// equal to = rho * [rho*i1]_t0d1 / rho
    /// equal to = i1t0d0 + rho * i1t0d1
    fn t0d1(&mut self, eta: f64) -> f64 {
        if eta != self.t0d1.0 {
            self.t0d1 = (
                eta,
                self.a0
                    + self.a1 * 2.0 * eta
                    + self.a2 * 3.0 * eta.powi(2)
                    + self.a3 * 4.0 * eta.powi(3)
                    + self.a4 * 5.0 * eta.powi(4)
                    + self.a5 * 6.0 * eta.powi(5)
                    + self.a6 * 7.0 * eta.powi(6),
            )
        }
        self.t0d1.1
    }
    /// equal to = rho^2 * [rho*i1]_t0d2 / rho
    /// equal to = 2.0 * rho * i1t0d1 + rho^2 * i1t0d2
    fn t0d2(&mut self, eta: f64) -> f64 {
        if eta != self.t0d2.0 {
            self.t0d2 = (
                eta,
                self.a1 * 2.0 * eta
                    + self.a2 * 6.0 * eta.powi(2)
                    + self.a3 * 12.0 * eta.powi(3)
                    + self.a4 * 20.0 * eta.powi(4)
                    + self.a5 * 30.0 * eta.powi(5)
                    + self.a6 * 42.0 * eta.powi(6),
            )
        }
        self.t0d2.1
    }
    /// equal to = rho^3 * [rho*i1]_t0d3 / rho
    /// equal to = 3.0 * rho^2 * i1t0d2 + rho^3 * i1t0d3
    fn t0d3(&mut self, eta: f64) -> f64 {
        if eta != self.t0d3.0 {
            self.t0d3 = (
                eta,
                self.a2 * 6.0 * eta.powi(2)
                    + self.a3 * 24.0 * eta.powi(3)
                    + self.a4 * 60.0 * eta.powi(4)
                    + self.a5 * 120.0 * eta.powi(5)
                    + self.a6 * 210.0 * eta.powi(6),
            )
        }
        self.t0d3.1
    }
    /// equal to = rho^4 * [rho*i1]_t0d4 / rho
    /// equal to = 4.0 * rho^3 * i1t0d3 + rho^4 * i1t0d4
    fn t0d4(&mut self, eta: f64) -> f64 {
        if eta != self.t0d4.0 {
            self.t0d4 = (
                eta,
                self.a3 * 24.0 * eta.powi(3)
                    + self.a4 * 120.0 * eta.powi(4)
                    + self.a5 * 360.0 * eta.powi(5)
                    + self.a6 * 840.0 * eta.powi(6),
            )
        }
        self.t0d4.1
    }
    /// equal to = { rho * [rho*i1]_t1d1 / rho } / eta1
    /// equal to = { i1t1d0 + rho * i1t1d1 } / eta1
    fn t1d1(&mut self, eta: f64) -> f64 {
        if eta != self.t1d1.0 {
            self.t1d1 = (
                eta,
                self.a1 * 2.0
                    + self.a2 * 6.0 * eta
                    + self.a3 * 12.0 * eta.powi(2)
                    + self.a4 * 20.0 * eta.powi(3)
                    + self.a5 * 30.0 * eta.powi(4)
                    + self.a6 * 42.0 * eta.powi(5),
            )
        }
        self.t1d1.1
    }
    /// equal to = { rho^2 * [rho*i1]_t1d2 / rho } / eta1
    /// equal to = { 2.0 * rho * i1t1d1 + rho^2 * i1t1d2 } / eta1
    fn t1d2(&mut self, eta: f64) -> f64 {
        if eta != self.t1d2.0 {
            self.t1d2 = (
                eta,
                self.a1 * 2.0
                    + self.a2 * 12.0 * eta
                    + self.a3 * 36.0 * eta.powi(2)
                    + self.a4 * 80.0 * eta.powi(3)
                    + self.a5 * 150.0 * eta.powi(4)
                    + self.a6 * 252.0 * eta.powi(5),
            )
        }
        self.t1d2.1
    }
    /// equal to = { rho^3 * [rho*i1]_t1d3 / rho } / eta1
    /// equal to = { 3.0 * rho^2 * i1t1d2 + rho^3 * i1t1d3 } / eta1
    fn t1d3(&mut self, eta: f64) -> f64 {
        if eta != self.t1d3.0 {
            self.t1d3 = (
                eta,
                self.a2 * 12.0 * eta
                    + self.a3 * 72.0 * eta.powi(2)
                    + self.a4 * 240.0 * eta.powi(3)
                    + self.a5 * 600.0 * eta.powi(4)
                    + self.a6 * 1260.0 * eta.powi(5),
            )
        }
        self.t1d3.1
    }
}
/// I2Term
struct I2Term {
    b0: f64,
    b1: f64,
    b2: f64,
    b3: f64,
    b4: f64,
    b5: f64,
    b6: f64,
    // cached variables 0
    eta0: (f64, f64),
    eta1: (f64, f64),
    eta2: (f64, f64),
    // cached variables 1
    t0d1: (f64, f64),
    t0d2: (f64, f64),
    t0d3: (f64, f64),
    t0d4: (f64, f64),
    t1d1: (f64, f64),
    t1d2: (f64, f64),
    t1d3: (f64, f64),
}
// Universal Model Constants
const B00: f64 = 0.7240946941;
const B01: f64 = 2.2382791861;
const B02: f64 = -4.0025849485;
const B03: f64 = -21.003576815;
const B04: f64 = 26.855641363;
const B05: f64 = 206.55133841;
const B06: f64 = -355.60235612;
const B10: f64 = -0.5755498075;
const B11: f64 = 0.6995095521;
const B12: f64 = 3.8925673390;
const B13: f64 = -17.215471648;
const B14: f64 = 192.67226447;
const B15: f64 = -161.82646165;
const B16: f64 = -165.20769346;
const B20: f64 = 0.0976883116;
const B21: f64 = -0.2557574982;
const B22: f64 = -9.1558561530;
const B23: f64 = 20.642075974;
const B24: f64 = -38.804430052;
const B25: f64 = 93.626774077;
const B26: f64 = -29.666905585;
impl I2Term {
    fn new(m: f64) -> Self {
        let m1 = (m - 1.0) / m; // (m-1)/m
        let m12 = (m - 1.0) * (m - 2.0) / m.powi(2); // (m-1)/m * (m-2)/m
        Self {
            b0: B00 + m1 * B10 + m12 * B20,
            b1: B01 + m1 * B11 + m12 * B21,
            b2: B02 + m1 * B12 + m12 * B22,
            b3: B03 + m1 * B13 + m12 * B23,
            b4: B04 + m1 * B14 + m12 * B24,
            b5: B05 + m1 * B15 + m12 * B25,
            b6: B06 + m1 * B16 + m12 * B26,
            // cached variables 0
            eta0: (0.0, 0.0),
            eta1: (0.0, 0.0),
            eta2: (0.0, 0.0),
            // cached variables 1
            t0d1: (0.0, 0.0),
            t0d2: (0.0, 0.0),
            t0d3: (0.0, 0.0),
            t0d4: (0.0, 0.0),
            t1d1: (0.0, 0.0),
            t1d2: (0.0, 0.0),
            t1d3: (0.0, 0.0),
        }
    }
    fn eta0(&mut self, eta: f64) -> f64 {
        if eta != self.eta0.0 {
            self.eta0 = (
                eta,
                self.b0
                    + self.b1 * eta
                    + self.b2 * eta.powi(2)
                    + self.b3 * eta.powi(3)
                    + self.b4 * eta.powi(4)
                    + self.b5 * eta.powi(5)
                    + self.b6 * eta.powi(6),
            )
        }
        self.eta0.1
    }
    fn eta1(&mut self, eta: f64) -> f64 {
        if eta != self.eta1.0 {
            self.eta1 = (
                eta,
                self.b1
                    + self.b2 * 2.0 * eta
                    + self.b3 * 3.0 * eta.powi(2)
                    + self.b4 * 4.0 * eta.powi(3)
                    + self.b5 * 5.0 * eta.powi(4)
                    + self.b6 * 6.0 * eta.powi(5),
            )
        }
        self.eta1.1
    }
    fn eta2(&mut self, eta: f64) -> f64 {
        if eta != self.eta2.0 {
            self.eta2 = (
                eta,
                self.b2 * 2.0
                    + self.b3 * 6.0 * eta
                    + self.b4 * 12.0 * eta.powi(2)
                    + self.b5 * 20.0 * eta.powi(3)
                    + self.b6 * 30.0 * eta.powi(4),
            )
        }
        self.eta2.1
    }
    /// equal to = [rho*i2]_t0d0 / rho
    /// equal to = i2t0d0
    /// equal to = eta0
    #[inline]
    fn t0d0(&mut self, eta: f64) -> f64 {
        self.eta0(eta)
    }
    /// equal to = { [rho*i2]_t1d0 / rho } / eta1
    /// equal to = { i2t1d0 } / eta1
    /// equal to = eta1
    #[inline]
    fn t1d0(&mut self, eta: f64) -> f64 {
        self.eta1(eta)
    }
    /// equal to = [rho*i2]_t2d0 / rho
    /// equal to = i2t2d0
    fn t2d0(&mut self, eta: f64, eta1: f64, eta2: f64) -> f64 {
        eta2 * self.eta1(eta) + eta1.powi(2) * self.eta2(eta)
    }
    /// equal to = rho * [rho*i2]_t0d1 / rho
    /// equal to = i2t0d0 + rho * i2t0d1
    fn t0d1(&mut self, eta: f64) -> f64 {
        if eta != self.t0d1.0 {
            self.t0d1 = (
                eta,
                self.b0
                    + self.b1 * 2.0 * eta
                    + self.b2 * 3.0 * eta.powi(2)
                    + self.b3 * 4.0 * eta.powi(3)
                    + self.b4 * 5.0 * eta.powi(4)
                    + self.b5 * 6.0 * eta.powi(5)
                    + self.b6 * 7.0 * eta.powi(6),
            )
        }
        self.t0d1.1
    }
    /// equal to = rho^2 * [rho*i2]_t0d2 / rho
    /// equal to = 2.0 * rho * i2t0d1 + rho^2 * i2t0d2
    fn t0d2(&mut self, eta: f64) -> f64 {
        if eta != self.t0d2.0 {
            self.t0d2 = (
                eta,
                self.b1 * 2.0 * eta
                    + self.b2 * 6.0 * eta.powi(2)
                    + self.b3 * 12.0 * eta.powi(3)
                    + self.b4 * 20.0 * eta.powi(4)
                    + self.b5 * 30.0 * eta.powi(5)
                    + self.b6 * 42.0 * eta.powi(6),
            )
        }
        self.t0d2.1
    }
    /// equal to = rho^3 * [rho*i2]_t0d3 / rho
    /// equal to = 3.0 * rho^2 * i2t0d2 + rho^3 * i2t0d3
    fn t0d3(&mut self, eta: f64) -> f64 {
        if eta != self.t0d3.0 {
            self.t0d3 = (
                eta,
                self.b2 * 6.0 * eta.powi(2)
                    + self.b3 * 24.0 * eta.powi(3)
                    + self.b4 * 60.0 * eta.powi(4)
                    + self.b5 * 120.0 * eta.powi(5)
                    + self.b6 * 210.0 * eta.powi(6),
            )
        }
        self.t0d3.1
    }
    /// equal to = rho^4 * [rho*i2]_t0d4 / rho
    /// equal to = 4.0 * rho^3 * i2t0d3 + rho^4 * i2t0d4
    fn t0d4(&mut self, eta: f64) -> f64 {
        if eta != self.t0d4.0 {
            self.t0d4 = (
                eta,
                self.b3 * 24.0 * eta.powi(3)
                    + self.b4 * 120.0 * eta.powi(4)
                    + self.b5 * 360.0 * eta.powi(5)
                    + self.b6 * 840.0 * eta.powi(6),
            )
        }
        self.t0d4.1
    }
    /// equal to = { rho * [rho*i2]_t1d1 / rho } / eta1
    /// equal to = { i2t1d0 + rho * i2t1d1 } / eta1
    fn t1d1(&mut self, eta: f64) -> f64 {
        if eta != self.t1d1.0 {
            self.t1d1 = (
                eta,
                self.b1 * 2.0
                    + self.b2 * 6.0 * eta
                    + self.b3 * 12.0 * eta.powi(2)
                    + self.b4 * 20.0 * eta.powi(3)
                    + self.b5 * 30.0 * eta.powi(4)
                    + self.b6 * 42.0 * eta.powi(5),
            )
        }
        self.t1d1.1
    }
    /// equal to = { rho^2 * [rho*i2]_t1d2 / rho } / eta1
    /// equal to = { 2.0 * rho * i2t1d1 + rho^2 * i2t1d2 } / eta1
    fn t1d2(&mut self, eta: f64) -> f64 {
        if eta != self.t1d2.0 {
            self.t1d2 = (
                eta,
                self.b1 * 2.0
                    + self.b2 * 12.0 * eta
                    + self.b3 * 36.0 * eta.powi(2)
                    + self.b4 * 80.0 * eta.powi(3)
                    + self.b5 * 150.0 * eta.powi(4)
                    + self.b6 * 252.0 * eta.powi(5),
            )
        }
        self.t1d2.1
    }
    /// equal to = { rho^3 * [rho*i2]_t1d3 / rho } / eta1
    /// equal to = { 3.0 * rho^2 * i2t1d2 + rho^3 * i2t1d3 } / eta1
    fn t1d3(&mut self, eta: f64) -> f64 {
        if eta != self.t1d3.0 {
            self.t1d3 = (
                eta,
                self.b2 * 12.0 * eta
                    + self.b3 * 72.0 * eta.powi(2)
                    + self.b4 * 240.0 * eta.powi(3)
                    + self.b5 * 600.0 * eta.powi(4)
                    + self.b6 * 1260.0 * eta.powi(5),
            )
        }
        self.t1d3.1
    }
}
