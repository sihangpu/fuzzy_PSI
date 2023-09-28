#![feature(test)]

extern crate test;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::RistrettoPoint;
use rand::rngs::OsRng;

fn fibonacci_u64(number: u64) -> u64 {
    let mut last: u64 = 1;
    let mut current: u64 = 0;
    let mut buffer: u64;
    let mut position: u64 = 1;

    return loop {
        if position == number {
            break current;
        }

        buffer = last;
        last = current;
        current = buffer + current;
        position += 1;
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[test]
    fn it_works() {
        assert_eq!(fibonacci_u64(1), 0);
        assert_eq!(fibonacci_u64(2), 1);
        assert_eq!(fibonacci_u64(12), 89);
        assert_eq!(fibonacci_u64(30), 514229);
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
