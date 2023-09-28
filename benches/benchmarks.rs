#![feature(test)]
extern crate test;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::RistrettoPoint;

use rand::rngs::OsRng;
use std::collections::HashMap;

extern crate f_psi;
use f_psi::okvs;

#[cfg(test)]
#[allow(unused_variables)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn okvs_works() {
        let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
        let n = 100;
        for j in 0..n {
            list.insert(j, RISTRETTO_BASEPOINT_POINT);
        }
        let m = (n as f64 * (1.0 + okvs::EPSILON)).ceil() as u64;
        let encoding = okvs::okvs_encode(&list);
        let decoding = okvs::okvs_decode(&encoding, 1);
        assert_eq!(encoding.len(), m as usize);
        assert_eq!(decoding, RISTRETTO_BASEPOINT_POINT);
    }

    #[bench]
    fn bench_pointadd(b: &mut Bencher) {
        let mut rng = OsRng;
        // Generate two random points on the curve.
        let point_ran1 = RistrettoPoint::random(&mut rng);
        let point_ran2 = RistrettoPoint::random(&mut rng);
        b.iter(|| {
            test::black_box(&point_ran1 + &point_ran2);
        });
    }
    #[bench]
    fn bench_scalarmult(b: &mut Bencher) {
        let scalar = Scalar::from(5u64);
        b.iter(|| {
            test::black_box(&scalar * &RISTRETTO_BASEPOINT_POINT);
        });
    }

    #[bench]
    fn bench_basemult(b: &mut Bencher) {
        let scalar = Scalar::from(5u64);
        b.iter(|| {
            test::black_box(&scalar * RISTRETTO_BASEPOINT_TABLE);
        });
    }
}
