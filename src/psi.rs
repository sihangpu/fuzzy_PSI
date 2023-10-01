use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use crate::okvs;
use fxhash::hash64;

pub const DIM: usize = 2;
pub const R: u64 = 30;
pub const BLK_CELLS: usize = 1 << DIM;
pub const SIDE_LEN: u64 = 2 * R;
pub type Point = [u64; DIM];

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
    _okvsgen: Vec<okvs::OkvsGen>,
}

impl Receiver {
    pub fn new(num_item: u64) -> Self {
        // generate random sk and pk
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng);
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE;
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(DIM);
        let mut _okvsgen: Vec<okvs::OkvsGen> = Vec::with_capacity(DIM);
        let window: usize = (2 * R + 1) as usize;
        let n: u64 = num_item * window as u64;
        for _ in 0..DIM {
            let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(n as usize);
            for _ in 0..n {
                let tem: Scalar = Scalar::random(&mut rng);
                pair.push((tem, tem * sk));
            }
            _pre_data.push(pair);
            _okvsgen.push(okvs::OkvsGen::new(n));
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

    pub fn get_windowsize(&self) -> usize {
        return self.window;
    }

    pub fn get_output_size_per_dim(&self) -> u64 {
        return self.n;
    }
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng();
        self._okvsgen.clear();
        self._pre_data.clear();
        for _ in 0..DIM {
            let mut pair = Vec::with_capacity(self.n as usize);
            for _ in 0..self.n {
                let tem = Scalar::random(&mut rng);
                pair.push((tem, tem * self.sk));
            }
            self._pre_data.push(pair);
            self._okvsgen.push(okvs::OkvsGen::new(self.n));
        }
    }
    pub fn msg(&mut self, pt_set: &Vec<Point>) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(DIM);
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        // for each dimension
        for i in 0..DIM {
            // for each receiver's point pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            {
                let blk = block(pt, SIDE_LEN, R);
                let key = hash64(&blk);

                // if pt_set[7][i] == pt[i] {
                //     println!("rec key at  7 : {:?}", blk);
                // }

                // for each possible value in [2R+1]
                let min = pt[i] - (self.window >> 1) as u64;
                for (j, pre_val) in pre_window.iter().enumerate() {
                    let key_ij = hash64(&(min + j as u64));
                    list.push((key ^ key_ij, *pre_val));
                }
            }
            result.push(self._okvsgen[i].encode(&list));
            list.clear();
        }
        self._okvsgen.clear();
        self._pre_data.clear();
        return result;
    }

    pub fn output(&self, msg_sender: &okvs::Encoding, window: usize) -> i32 {
        let mut count = 0;
        for values in msg_sender.windows(window).step_by(window) {
            for (u, v) in values.iter() {
                if (self.sk * u) == *v {
                    count += 1;
                    break;
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

    pub fn get_output_size(&self) -> u64 {
        return self.m;
    }

    pub fn get_windowsize(&self) -> usize {
        return self.window;
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
    pub fn msg(&mut self, encodings: &Vec<okvs::Encoding>, pt_set: &Vec<Point>) -> okvs::Encoding {
        let mut blk: Point = [0u64; DIM];
        let mut uv: okvs::Encoding =
            vec![(RistrettoPoint::identity(), RistrettoPoint::identity()); self.m as usize];
        let mut tem: (RistrettoPoint, RistrettoPoint);
        // for each senders's point pt
        for (ind, (pt, coin_window)) in pt_set
            .iter()
            .zip(self._coins.windows(self.window).step_by(self.window))
            .enumerate()
        {
            let cel = cell(pt, SIDE_LEN);
            // for each possible block
            for (i, coins) in coin_window.iter().enumerate() {
                let uv_i = ind * self.window + i;
                for j in 0..DIM {
                    if (i >> j) & 1 == 1 {
                        blk[j] = cel[j] - 1;
                    } else {
                        blk[j] = cel[j];
                    }
                }
                // if ind == 9 {
                //     println!("key at {} in {}-blk: {:?}", ind, i, blk);
                // }
                let key = hash64(&blk);
                // for each dimension
                for j in 0..DIM {
                    let key_ij = hash64(&(pt[j] as u64));
                    tem = okvs::okvs_decode(&encodings[j], key ^ key_ij);
                    uv[uv_i].0 += tem.0;
                    uv[uv_i].1 += tem.1;
                }
                // finalize
                uv[uv_i].0 = coins.2 * uv[uv_i].0 + coins.0;
                uv[uv_i].1 = coins.2 * uv[uv_i].1 + coins.1;
            }
        }
        self._coins.clear();
        return uv;
    }
}
