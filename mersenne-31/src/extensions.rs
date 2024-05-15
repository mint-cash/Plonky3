use core::mem;
use p3_field::{
    AbstractField, BinomialExtensionAlgebra, BinomialExtensionParams, Complex, ComplexExtendable,
};

use crate::Mersenne31;

const fn m31(x: u32) -> Mersenne31 {
    Mersenne31::new(x)
}
const fn m31s<const N: usize>(xs: [u32; N]) -> [Mersenne31; N] {
    let mut ys = [Mersenne31::new(0); N];
    let mut i = 0;
    while i < N {
        ys[i] = m31(xs[i]);
        i += 1;
    }
    ys
}

impl ComplexExtendable for Mersenne31 {
    const COMPLEX_GEN: [Self; 2] = m31s([12, 1]);
    const CIRCLE_TWO_ADICITY: usize = 31;
    fn circle_two_adic_generator(bits: usize) -> [Self; 2] {
        let base = Complex::new(m31(311_014_874), m31(1_584_694_829));
        base.exp_power_of_2(Self::CIRCLE_TWO_ADICITY - bits).0
    }
}

#[cfg(test)]
mod test_m31_complex {
    use super::*;
    use p3_field::Complex;
    use p3_field_testing::test_field;
    test_field!(Complex<Mersenne31>);
}

pub type Mersenne31Cubic = BinomialExtensionAlgebra<Mersenne31, 3, Mersenne31CubicParams>;

#[derive(Debug)]
pub struct Mersenne31CubicParams;

impl BinomialExtensionParams<Mersenne31, 3> for Mersenne31CubicParams {
    const W: Mersenne31 = m31(5);
    const ORDER_D_SUBGROUP: [Mersenne31; 3] = m31s([1, 1513477735, 634005911]);
    const GEN: [Mersenne31; 3] = m31s([10, 1, 0]);
}

#[cfg(test)]
mod test_m31_cubic {
    use super::*;
    use p3_field::Extension;
    use p3_field_testing::test_field;
    test_field!(Extension<Mersenne31Cubic>);
}

/*
impl BinomiallyExtendable<3> for Mersenne31 {
    // ```sage
    // p = 2^31 - 1
    // F = GF(p)
    // R.<x> = F[]
    // assert (x^3 - 5).is_irreducible()
    // ```
    fn w() -> Self {
        Self::new(5)
    }
    // ```sage
    // F(5)^((p-1)/3)
    // ```
    fn dth_root() -> Self {
        Self::new(1513477735)
    }
    // ```sage
    // F.extension(x^3 - 5, 'u').multiplicative_generator()
    // ```
    fn ext_generator() -> [Self; 3] {
        [Self::new(10), Self::new(1), Self::zero()]
    }
}

impl HasComplexBinomialExtension<2> for Mersenne31 {
    // Verifiable in Sage with
    // ```sage
    // p = 2**31 - 1  # Mersenne31
    // F = GF(p)  # The base field GF(p)
    // R.<x> = F[]  # The polynomial ring over F
    // K.<i> = F.extension(x^2 + 1)  # The complex extension field
    // R2.<y> = K[]
    // f2 = y^2 - i - 2
    // assert f2.is_irreducible()
    // ```
    fn w() -> Complex<Self> {
        Complex::new(Mersenne31::new(2), Mersenne31::one())
    }

    // Verifiable in Sage with
    // ```sage
    // K2.<j> = K.extension(f2)
    //  g = j + 6
    // for f in factor(p^4 - 1):
    //   assert g^((p^4-1) // f) != 1
    // ```
    fn ext_generator() -> [Complex<Self>; 2] {
        [Complex::new_real(Mersenne31::new(6)), Complex::one()]
    }

    // DTH_ROOT = W^((p^2 - 1)/2).
    fn dth_root() -> Complex<Self> {
        Complex::new_real(Mersenne31::new(2147483646))
    }
}

impl HasTwoAdicComplexBinomialExtension<2> for Mersenne31 {
    const COMPLEX_EXT_TWO_ADICITY: usize = 33;

    fn complex_ext_two_adic_generator(bits: usize) -> [Complex<Self>; 2] {
        assert!(bits <= 33);
        if bits == 33 {
            [
                Complex::zero(),
                Complex::new(Mersenne31::new(1437746044), Mersenne31::new(946469285)),
            ]
        } else {
            [Complex::two_adic_generator(bits), Complex::zero()]
        }
    }
}

impl HasComplexBinomialExtension<3> for Mersenne31 {
    // Verifiable in Sage with
    // ```sage
    // p = 2**31 - 1  # Mersenne31
    // F = GF(p)  # The base field GF(p)
    // R.<x> = F[]  # The polynomial ring over F
    // K.<i> = F.extension(x^2 + 1)  # The complex extension field
    // R2.<y> = K[]
    // f2 = y^3 - 5*i
    // assert f2.is_irreducible()
    // ```
    fn w() -> Complex<Self> {
        Complex::new_imag(Mersenne31::new(5))
    }

    // DTH_ROOT = W^((p^2 - 1)/2).
    fn dth_root() -> Complex<Self> {
        Complex::new_real(Mersenne31::new(634005911))
    }

    // Verifiable in Sage with
    // ```sage
    // K2.<j> = K.extension(f2)
    //  g = j + 5
    // for f in factor(p^6 - 1):
    //   assert g^((p^6-1) // f) != 1
    // ```
    fn ext_generator() -> [Complex<Self>; 3] {
        [
            Complex::new_real(Mersenne31::new(5)),
            Complex::new_real(Mersenne31::one()),
            Complex::zero(),
        ]
    }
}

impl HasTwoAdicComplexBinomialExtension<3> for Mersenne31 {
    const COMPLEX_EXT_TWO_ADICITY: usize = 32;

    fn complex_ext_two_adic_generator(bits: usize) -> [Complex<Self>; 3] {
        field_to_array::<Complex<Self>, 3>(Complex::two_adic_generator(bits))
    }
}

#[cfg(test)]
mod test_cubic_extension {
    use p3_field::extension::{BinomialExtensionField, Complex};
    use p3_field_testing::{test_field, test_two_adic_extension_field};

    use crate::Mersenne31;

    type F = Complex<Mersenne31>;
    type EF = BinomialExtensionField<F, 3>;

    test_field!(super::EF);

    test_two_adic_extension_field!(super::F, super::EF);
}

#[cfg(test)]
mod test_quadratic_extension {

    use p3_field::extension::{BinomialExtensionField, Complex};
    use p3_field_testing::{test_field, test_two_adic_extension_field};

    use crate::Mersenne31;

    type F = Complex<Mersenne31>;
    type EF = BinomialExtensionField<F, 2>;

    test_field!(super::EF);

    test_two_adic_extension_field!(super::F, super::EF);
}

*/