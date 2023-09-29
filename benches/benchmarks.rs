#![feature(test)]
#![allow(unused_imports)]
#![allow(dead_code)]

extern crate test;

extern crate f_psi;
use f_psi::okvs;

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::time::Instant;
    use test::Bencher;

    use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
    use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
    use curve25519_dalek::scalar::Scalar;
    use curve25519_dalek::RistrettoPoint;

    use fxhash::hash64;
    use rand::rngs::OsRng;
    use std::collections::HashMap;

    // #[test]
    // fn hash_test() {
    //     //test hash64
    //     let h1 = hash64(&(0 ^ 0x1234567890abcdefu64));
    //     let h1c = hash64(&(0 ^ 0x1234567890abcdefu64));
    //     let seed = hash64(&(0 ^ 0x1234567890abcdefu64));
    //     let hash_band: Vec<Scalar> = (0..40).map(|i| Scalar::from((seed >> i) & 0x01)).collect();
    //     //print hash_band
    //     for i in 0..40 {
    //         print!("{:?} ", hash_band[i] == Scalar::ONE);
    //     }
    //     println!("h1: {:x}", h1);
    //     assert_eq!(h1, h1c);
    // }
    //
    // #[test]
    // fn gbf_test() {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let N = 100;
    //     for j in 0..N {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(N);
    //     gbf.encode(&list);
    //     let decoding = gbf.decode(12);
    //     assert_eq!(decoding, RISTRETTO_BASEPOINT_POINT);
    // }

    const N: u64 = 1000;

    #[test]
    fn okvs_test() {
        let mut list: HashMap<u64, (Scalar, Scalar)> = HashMap::new();
        for j in 0..N {
            list.insert(j + 1, (Scalar::ONE, Scalar::ONE));
        }
        let mut okvsmod = okvs::OKVS::new(N);

        let now = Instant::now();

        test::black_box(okvsmod.encode(&list));

        let elapsed = now.elapsed();
        println!("{} items, Elapsed Enc: {:.2?}", N, elapsed);
        let decoding = okvsmod.decode(39);
        assert_eq!(
            decoding,
            (RISTRETTO_BASEPOINT_POINT, RISTRETTO_BASEPOINT_POINT)
        );
    }

    #[bench]
    fn bench_okvs_encode(b: &mut Bencher) {
        let mut list: HashMap<u64, (Scalar, Scalar)> = HashMap::new();
        println!("{} items, OKVS.Encode", N);
        for j in 0..N {
            list.insert(j + 1, (Scalar::ONE, Scalar::ONE));
        }
        let mut okvsmod = okvs::OKVS::new(N);
        b.iter(|| {
            test::black_box(okvsmod.encode(&list));
        });
    }
    #[bench]
    fn bench_okvs_decode(b: &mut Bencher) {
        let mut list: HashMap<u64, (Scalar, Scalar)> = HashMap::new();
        println!("{} items, OKVS.Decode", N);
        for j in 0..N {
            list.insert(j + 1, (Scalar::ONE, Scalar::ONE));
        }
        let mut okvsmod = okvs::OKVS::new(N);
        okvsmod.encode(&list);
        b.iter(|| {
            test::black_box(okvsmod.decode(39));
        });
    }

    // #[bench]
    // fn bench_gbf_encode(b: &mut Bencher) {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let N = 2000;
    //     for j in 0..N {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(N);
    //     b.iter(|| {
    //         test::black_box(gbf.encode(&list));
    //     });
    // }
    // #[bench]
    // fn bench_gbf_decode(b: &mut Bencher) {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let N = 2000;
    //     for j in 0..N {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(N);
    //     b.iter(|| {
    //         test::black_box(gbf.decode(8));
    //     });
    // }

    // fn add_n(p1: &RistrettoPoint, p2: &RistrettoPoint, N: usize) -> RistrettoPoint {
    //     let mut p: RistrettoPoint = RISTRETTO_BASEPOINT_POINT;
    //     for _ in 0..N {
    //         test::black_box(p = p1 - p2);
    //     }
    //     hash64(&9238u64);
    //     return p;
    // }
    // fn vec_it(v: &Vec<RistrettoPoint>) {
    //     for i in 0..v.len() {
    //         test::black_box(v[i]);
    //     }
    // }
    // #[bench]
    // fn bench_pointadd(b: &mut Bencher) {
    //     let mut rng = OsRng;
    //     // Generate two random points on the curve.
    //     let point_ran1 = RistrettoPoint::random(&mut rng);
    //     let point_ran2 = RistrettoPoint::random(&mut rng);
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();

    //     let N = 2000;
    //     let m = (N as f64 * okvs::FACTOR).floor() as u64;
    //     for j in 0..N {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut enc: Vec<RistrettoPoint> =
    //         (0..m).map(|_| RistrettoPoint::random(&mut OsRng)).collect();
    //     b.iter(|| {
    //         // test::black_box(add_n(&point_ran1, &point_ran2, 2000 * 40));
    //         vec_it(&enc);
    //         // test::black_box(okvs::gbf_encode(&list, &mut enc));
    //     });
    // }

    // fn scalar_fp(s: &Scalar, t: &Scalar, N: usize) {
    //     let mut r = Scalar::ZERO;
    //     for _ in 0..N {
    //         r = s.invert();
    //         for _ in 0..40 {
    //             r += r * t;
    //         }
    //     }
    // }
    // #[bench]
    // fn bench_scalar(b: &mut Bencher) {
    //     let scalar1 = Scalar::from(5u64);
    //     let scalar2 = Scalar::from(7u64);
    //     b.iter(|| {
    //         // test::black_box(&scalar1.invert());
    //         test::black_box(scalar1 == Scalar::ZERO);
    //         // test::black_box(scalar_fp(&scalar1, &scalar2, 2000));
    //     });
    // }

    // #[bench]
    // fn bench_scalarmult(b: &mut Bencher) {
    //     let scalar = Scalar::from(5u64);
    //     b.iter(|| {
    //         test::black_box(&scalar * &RISTRETTO_BASEPOINT_POINT);
    //     });
    // }

    // #[bench]
    // fn bench_basemult(b: &mut Bencher) {
    //     let scalar = Scalar::from(5u64);
    //     b.iter(|| {
    //         test::black_box(&scalar * RISTRETTO_BASEPOINT_TABLE);
    //     });
    // }
}
