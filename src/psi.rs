use std::collections::HashMap;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use crate::okvs::OKVS;
use fxhash::hash64;
const DIM: usize = 2;
const R: u64 = 20;

const BLK_CELLS: usize = 2usize.pow(DIM as u32);
const SHIFT_ORIGIN: u64 = 2 * R;

type Point = [u64; DIM];

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

struct Receiver {
    n: u64,
    pk: RistrettoPoint,
    sk: Scalar,
    _pre_data: Vec<(Scalar, Scalar)>,
    _okvsgen: OKVS,
}

impl Receiver {
    pub fn new(num_item: u64) -> Receiver {
        // generate random sk and pk
        let mut rng = rand::thread_rng();
        let sk = Scalar::random(&mut rng);
        let pk = &sk * RISTRETTO_BASEPOINT_TABLE;
        let mut _pre_data: Vec<(Scalar, Scalar)> = Vec::with_capacity(num_item as usize);
        for _ in 0..num_item {
            let tem = Scalar::random(&mut rng);
            _pre_data.push((tem, tem * sk));
        }
        return Receiver {
            n: num_item,
            pk,
            sk,
            _pre_data,
            _okvsgen: OKVS::new(num_item),
        };
    }
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng();
        self._pre_data.clear();
        for _ in 0..self.n {
            let tem = Scalar::random(&mut rng);
            self._pre_data.push((tem, tem * self.sk));
        }
    }
    pub fn receive(
        &mut self,
        pt_set: Vec<Point>,
        radius: u64,
    ) -> Vec<(RistrettoPoint, RistrettoPoint)> {
        let mut list: HashMap<u64, (Scalar, Scalar)> = HashMap::new();
        for (i, pt) in pt_set.iter().enumerate() {
            let blk = block(&pt, 2 * R, R);
            let key = hash64(&blk);
            list.insert(key, self._pre_data[i]);
        }
        return self._okvsgen.encode(&list);
    }
}

#[test]
fn main() {
    println!("2R: 40,  cell: {:?}", cell(&[40, 80], 2 * R));
    println!("2R: 40,  cell: {:?}", block(&[50, 40], 2 * R, R));
}
