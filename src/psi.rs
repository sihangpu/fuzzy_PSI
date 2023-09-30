use std::collections::HashMap;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;
use rand::Rng;

use crate::okvs::okvs_decode;
use crate::okvs::OKVS;
use fxhash::hash64;

const DIM: usize = 2;
const R: u64 = 5;

const BLK_CELLS: usize = 1 << DIM;
const SHIFT_ORIGIN: u64 = 2 * R;

type Point = [u64; DIM];
type Encoding = Vec<(RistrettoPoint, RistrettoPoint)>;

fn sample_test_data_points(num: usize) -> Vec<Point> {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Point> = Vec::with_capacity(num);
    for _ in 0..num {
        let mut point: Point = [0u64; DIM];
        for i in 0..DIM {
            point[i] = rng.gen_range(SHIFT_ORIGIN..=(1 << 31));
        }
        points.push(point);
    }
    return points;
}

#[inline]
fn cell(p: &Point, sidele: u64) -> Point {
    let mut bot_left_corner: Point = [0u64; DIM];
    for i in 0..DIM {
        bot_left_corner[i] = p[i] / sidele;
    }
    return bot_left_corner;
}

#[inline]
fn block(p: &Point, sidele: u64, radius: u64) -> Point {
    let mut min: Point = [0u64; DIM];
    for i in 0..DIM {
        min[i] = p[i] - radius;
    }
    return cell(&min, sidele);
}

pub struct Receiver {
    window: usize,
    n: u64,
    pk: RistrettoPoint,
    sk: Scalar,
    _pre_data: Vec<Vec<(Scalar, Scalar)>>,
    _okvsgen: Vec<OKVS>,
}

impl Receiver {
    pub fn new(num_item: u64) -> Self {
        // generate random sk and pk
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng);
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE;
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(DIM);
        let mut _okvsgen: Vec<OKVS> = Vec::with_capacity(DIM);
        let window: usize = (2 * R + 1) as usize;
        let n: u64 = num_item * window as u64;
        for _ in 0..DIM {
            let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(n as usize);
            for _ in 0..n {
                let tem: Scalar = Scalar::random(&mut rng);
                pair.push((tem, tem * sk));
            }
            _pre_data.push(pair);
            _okvsgen.push(OKVS::new(n));
        }
        return Receiver {
            window,
            n,
            pk,
            sk,
            _pre_data,
            _okvsgen,
        };
    }

    pub fn publish_pk(&self) -> RistrettoPoint {
        return self.pk;
    }
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng();
        self._okvsgen.clear();
        self._pre_data.clear();
        for i in 0..DIM {
            let mut pair = Vec::with_capacity(self.n as usize);
            for _ in 0..self.n {
                let tem = Scalar::random(&mut rng);
                pair.push((tem, tem * self.sk));
            }
            self._pre_data.push(pair);
            self._okvsgen.push(OKVS::new(self.n));
        }
    }
    pub fn msg(&mut self, pt_set: &Vec<Point>) -> Vec<Encoding> {
        let mut result: Vec<Encoding> = Vec::with_capacity(DIM);
        // for each dimension
        for i in 0..DIM {
            let mut list: HashMap<u64, (Scalar, Scalar)> = HashMap::new();
            // for each receiver's point pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            {
                let blk = block(pt, R << 1, R);
                let key = hash64(&blk);
                // for each possible value in [2R+1]
                let min = pt[i] - (self.window >> 1) as u64;
                for (j, pre_val) in pre_window.iter().enumerate() {
                    let key_ij = hash64(&(min + j as u64));
                    list.insert(key ^ key_ij, *pre_val);
                }
            }
            result.push(self._okvsgen[i].encode(&list));
        }
        self._okvsgen.clear();
        self._pre_data.clear();
        return result;
    }

    pub fn output(&self, msg_sender: &Encoding, window: usize) -> i32 {
        let mut count = 0;
        for values in msg_sender.windows(window).step_by(window) {
            for (u, v) in values.iter() {
                if self.sk * u == *v {
                    count += 1;
                }
            }
        }
        return count;
    }
}

pub struct Sender {
    m: u64,
    window: usize,
    pk: RistrettoPoint,
    _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)>,
}
impl Sender {
    pub fn new(num_item: u64, pk_rec: RistrettoPoint) -> Self {
        // generate random sk and pk
        let window = BLK_CELLS;
        let m = num_item * window as u64;
        let mut rng = rand::thread_rng();
        let pk = pk_rec;
        let mut _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)> =
            Vec::with_capacity(m as usize);
        for _ in 0..m {
            let a = Scalar::random(&mut rng);
            _coins.push((
                &a * RISTRETTO_BASEPOINT_TABLE,
                a * pk_rec,
                Scalar::random(&mut rng),
            ));
        }
        return Sender {
            m,
            window,
            pk,
            _coins,
        };
    }
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng();
        self._coins.clear();
        for _ in 0..self.m {
            let a = Scalar::random(&mut rng);
            self._coins.push((
                &a * RISTRETTO_BASEPOINT_TABLE,
                a * self.pk,
                Scalar::random(&mut rng),
            ));
        }
    }
    pub fn msg(&mut self, encodings: &Vec<Encoding>, pt_set: &Vec<Point>) -> Encoding {
        let mut blk: Point = [0u64; DIM];
        let mut uv: Encoding =
            vec![(RistrettoPoint::identity(), RistrettoPoint::identity()); self.m as usize];
        let mut tem: (RistrettoPoint, RistrettoPoint);
        // for each senders's point pt
        for (ind, (pt, coin_window)) in pt_set
            .iter()
            .zip(self._coins.windows(self.window).step_by(self.window))
            .enumerate()
        {
            let cel = cell(pt, R << 1);
            // for each possible block
            for (i, coins) in coin_window.iter().enumerate() {
                for j in 0..DIM {
                    if (i >> j) & 1 == 1 {
                        blk[j] = cel[j] - R;
                    } else {
                        blk[j] = cel[j];
                    }
                }
                let key = hash64(&blk);
                // for each dimension
                for j in 0..DIM {
                    let key_ij = hash64(&(pt[j] as u64));
                    tem = okvs_decode(&encodings[j], key ^ key_ij);
                    uv[ind].0 += tem.0;
                    uv[ind].1 += tem.1;
                }
                // finalize
                uv[ind].0 = coins.2 * uv[ind].0 + coins.0;
                uv[ind].1 = coins.2 * uv[ind].1 + coins.1;
            }
        }
        self._coins.clear();
        return uv;
    }
}
#[test]
fn main() {
    println!("2R: 40,  cell: {:?}", cell(&[40, 80], 2 * R));
    println!("2R: 40,  cell: {:?}", block(&[50, 40], 2 * R, R));
    println!("{:?}", sample_test_data_points(10));
}
