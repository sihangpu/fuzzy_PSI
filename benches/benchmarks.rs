#![feature(test)]
#![allow(unused_imports)]
#![allow(dead_code)]

extern crate f_psi;
use f_psi::okvs;
use f_psi::protocol;
use f_psi::psi;

#[cfg(test)]
extern crate test;
mod tests {
    use super::*;
    use rand::rngs::OsRng;
    use rand::Rng;
    use std::time::Instant;
    use std::vec;
    use test::Bencher;

    use std::collections::HashSet;

    use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
    use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
    use curve25519_dalek::scalar::Scalar;
    use curve25519_dalek::traits::Identity;
    use curve25519_dalek::RistrettoPoint;

    use fxhash::hash64;

    fn sample_test_data_points(num: usize) -> Vec<psi::Point> {
        let mut rng = rand::thread_rng();
        let mut points: Vec<psi::Point> = Vec::with_capacity(num);
        for _ in 0..num {
            let mut point: psi::Point = [0u64; psi::DIM];
            for i in 0..psi::DIM {
                point[i] = rng.gen_range(psi::SIDE_LEN..=(1 << 31));
            }
            points.push(point);
        }
        return points;
    }

    #[test]
    fn protocol() {
        let n = 2000;
        let m = 10000;
        let data_r = sample_test_data_points(n);
        let mut data_s = sample_test_data_points(m);
        data_s[9][0] = data_r[7][0] - 15;
        data_s[9][1] = data_r[7][1] + 15;
        data_s[11][0] = data_r[7][0] + psi::R;
        data_s[11][1] = data_r[7][1] - psi::R;
        let (rcr, sdr) = protocol::setup(n, m, false, 0);
        let now = Instant::now();
        test::black_box(protocol::run_standard(rcr, sdr, data_r, data_s));
        let elapsed = now.elapsed();
        println!("Elapsed Time for Protocol: {:.2?}", elapsed);
    }

    #[test]
    fn protocol_apart() {
        let n = 2000;
        let m = 1000;
        let data_r = sample_test_data_points(n);
        let mut data_s = sample_test_data_points(m);
        data_s[9][0] = data_r[7][0] - 15;
        data_s[9][1] = data_r[7][1] + 15;
        data_s[11][0] = data_r[7][0] + psi::R;
        data_s[11][1] = data_r[7][1] - psi::R;
        let (rcr, sdr) = protocol::setup(n, m, true, 0);
        let now = Instant::now();
        test::black_box(protocol::run_standard_apart(rcr, sdr, data_r, data_s));
        let elapsed = now.elapsed();
        println!("Elapsed Time for Protocol: {:.2?}", elapsed);
    }

    #[test]
    fn protocol_lp() {
        let n = 2000;
        let m = 1000;
        let data_r = sample_test_data_points(n);
        let mut data_s = sample_test_data_points(m);
        data_s[9][0] = data_r[7][0] - 15;
        data_s[9][1] = data_r[7][1] + 15;
        data_s[11][0] = data_r[7][0] + psi::R;
        data_s[11][1] = data_r[7][1] - psi::R;
        let (rcr, sdr) = protocol::setup(n, m, true, 2);
        let now = Instant::now();
        test::black_box(protocol::run_standard_lp(rcr, sdr, data_r, data_s, 2));
        let elapsed = now.elapsed();
        println!("Elapsed Time for Protocol: {:.2?}", elapsed);
    }

    #[test]
    fn psi_lp() {
        let n = 2000;
        let m = 1000000;
        let data_r = sample_test_data_points(n);
        let mut data_s = sample_test_data_points(m);
        data_s[9][0] = data_r[7][0] - 15;
        data_s[9][1] = data_r[7][1] + 15;
        data_s[11][0] = data_r[7][0] + psi::R;
        data_s[11][1] = data_r[7][1] - psi::R;
        // println!("data_r points: {:?}", data_r[7]);
        // println!("data_s points: {:?}", data_s[9]);

        let mut rec_instance = psi::Receiver::new(n as u64, true);
        let send_instance = psi::Sender::new(m as u64, rec_instance.publish_pk(), true, 2);

        let now = Instant::now();
        let msg1 = rec_instance.lp_msg_apart(&data_r, 2);

        let elapsed = now.elapsed();
        println!(
            "{} items, Elapsed Time for Encoding (optimize=0): {:.2?}",
            rec_instance.get_output_size_per_dim() * psi::DIM as u64,
            elapsed
        );

        let sendnow = Instant::now();
        let mut msgvec: Vec<(okvs::PointPair, HashSet<u32>)> = Vec::with_capacity(m);
        for i in 0..m {
            msgvec.push(send_instance.lp_send_msg_single_apart(&msg1, &data_s[i], i));
        }
        let sendelapsed = sendnow.elapsed();
        println!(
            "{} items, Elapsed Time for Decoding (optimize=0): {:.2?}",
            send_instance.get_output_size(),
            sendelapsed
        );

        let mut c = 0u32;

        let recnow = Instant::now();
        // let out = rec_instance.output(&msg2, send_instance.get_windowsize());
        for i in 0..m {
            c += rec_instance.lp_post_process(&msgvec[i]);
        }
        let recoutputelapsed = recnow.elapsed();
        println!(
            "{} items, Elapsed Time for Finishing (optimize=0): {:.2?}",
            send_instance.get_output_size(),
            recoutputelapsed
        );
        println!("Total : {:.2?}", now.elapsed());
        println!("out: {}", c);
    }

    #[test]
    fn okvs_test() {
        let n = 2000 * 480;
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        for j in 0..n {
            list.push((
                hash64(&j),
                (Scalar::ONE, Scalar::from(100u64) * Scalar::ONE),
            ));
        }

        let mut okvs_instance = okvs::OkvsGen::new(n);

        let now = Instant::now();

        let data = okvs_instance.encode(&list);

        let elapsed = now.elapsed();
        println!(
            "{} items, Elapsed Time for Encoding (optimize=0): {:.2?}",
            n, elapsed
        );

        for i in 0..n {
            let decoding = okvs::okvs_decode(&data, hash64(&i));
            assert_eq!(
                decoding,
                (
                    RISTRETTO_BASEPOINT_POINT,
                    Scalar::from(100u64) * RISTRETTO_BASEPOINT_POINT
                )
            );
        }
    }

    #[bench]
    fn bench_okvs_encode(b: &mut Bencher) {
        let n = 2000;
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        println!("{} items, OkvsGen.Encode", n);
        for j in 0..n {
            list.push((hash64(&j), (Scalar::ONE, Scalar::ONE)));
        }
        let mut okvsmod = okvs::OkvsGen::new(n);
        b.iter(|| {
            test::black_box(okvsmod.encode(&list));
        });
    }
    #[inline]
    fn dec(dat: &okvs::Encoding, len: u64) {
        for i in 0..len {
            test::black_box(okvs::okvs_decode(dat, i));
        }
    }
    #[bench]
    fn bench_okvs_decode(b: &mut Bencher) {
        let n = 2000;
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        println!("{} items, OkvsGen.Decode", n);
        for j in 0..n {
            list.push((hash64(&j), (Scalar::ONE, Scalar::ONE)));
        }
        let mut okvsmod = okvs::OkvsGen::new(n);
        let data = okvsmod.encode(&list);
        let keys: Vec<u64> = (0..n).collect();
        b.iter(|| okvs::okvs_decode_batch(&data, &keys));
    }

    //----------------- gbf test -----------------
    // #[test]
    // fn gbf_test() {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let n = 100;
    //     for j in 0..n {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(n);
    //     gbf.encode(&list);
    //     let decoding = gbf.decode(12);
    //     assert_eq!(decoding, RISTRETTO_BASEPOINT_POINT);
    // }

    // #[bench]
    // fn bench_gbf_encode(b: &mut Bencher) {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let n = 2000;
    //     for j in 0..n {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(n);
    //     b.iter(|| {
    //         test::black_box(gbf.encode(&list));
    //     });
    // }
    // #[bench]
    // fn bench_gbf_decode(b: &mut Bencher) {
    //     let mut list: HashMap<u64, RistrettoPoint> = HashMap::new();
    //     let n = 2000;
    //     for j in 0..n {
    //         list.insert(j, RISTRETTO_BASEPOINT_POINT);
    //     }
    //     let mut gbf = okvs::GBF::new(n);
    //     b.iter(|| {
    //         test::black_box(gbf.decode(8));
    //     });
    // }
    //
    //----------------- elliptic_curve test -----------------
    // fn add_n(
    //     p1: &RistrettoPoint,
    //     p2: &RistrettoPoint,
    //     n: usize,
    //     v: &okvs::Encoding,
    //     ind: &Vec<usize>,
    // ) -> RistrettoPoint {
    //     let mut p: RistrettoPoint = RISTRETTO_BASEPOINT_POINT;
    //     for i in 0..n {
    //         test::black_box(p = v[ind[i]].0 - p1);
    //         test::black_box(p = v[ind[i]].1 - p2);
    //     }
    //     return p;
    // }

    // #[bench]
    // fn bench_pointadd(b: &mut Bencher) {
    //     let mut rng = OsRng;
    //     let n: usize = 2000*40;
    //     // Generate two random points on the curve.
    //     let point_ran1 = RistrettoPoint::random(&mut rng);
    //     let point_ran2 = RistrettoPoint::random(&mut rng);
    //     let vec: okvs::Encoding = (0..n)
    //         .map(|_| {
    //             (
    //                 RistrettoPoint::random(&mut rng),
    //                 RistrettoPoint::random(&mut rng),
    //             )
    //         })
    //         .collect();
    //     let index: Vec<usize> = (0..n).map(|_| rng.gen_range(0..n)).collect();
    //     b.iter(|| {
    //         test::black_box(add_n(&point_ran1, &point_ran2, n, &vec, &index));
    //         // vec_it(&enc);
    //         // test::black_box(okvs::gbf_encode(&list, &mut enc));
    //     });
    // }

    // fn scalar_fp(s: &Scalar, t: &Scalar) {
    //     let mut r = Scalar::ZERO;
    //     for _ in 0..n {
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
    //         // test::black_box(scalar1 == Scalar::ZERO);
    //         test::black_box(scalar_fp(&scalar1, &scalar2));
    //     });
    // }

    // #[bench]
    // fn bench_scalarmult(b: &mut Bencher) {
    //     let scalar = Scalar::from(0u64);
    //     b.iter(|| {
    //         test::black_box(assert_eq!(
    //             scalar * &RISTRETTO_BASEPOINT_POINT,
    //             RistrettoPoint::identity()
    //         ));
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

//----------------- criterion -----------------
// use criterion::{black_box, criterion_group, criterion_main, Criterion};
// fn bench_okvs_encode(b: &mut Criterion) {
//     let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
//     println!("{} items, OkvsGen.Encode", n);
//     for j in 0..n {
//         list.push((j + 1, (Scalar::ONE, Scalar::ONE)));
//     }
//     let mut okvs_instance = okvs::OkvsGen::new(n);
//     b.bench_function("okvs_encode", |b| {
//         b.iter(|| {
//             black_box(okvs_instance.encode(&list));
//         })
//     });
// }
// criterion_group!(benches, bench_okvs_encode);
// criterion_main!(benches);
//----------------- criterion -----------------
