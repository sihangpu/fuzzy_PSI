use std::collections::HashSet;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use crate::okvs;
use fxhash::hash64;

pub const DIM: usize = 2;
pub const R: u64 = 20;
pub const BLK_CELLS: usize = 1 << DIM;
pub const SIDE_LEN: u64 = 2 * R;
pub const R_L2: u64 = R * R;
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

#[inline]
fn l2_dist(p1: &Point, p2: &Point) -> u64 {
    let mut sum: u64 = 0;
    let mut diff: u64;
    for i in 0..DIM {
        diff = (p1[i] as i64 - p2[i] as i64).abs() as u64;
        sum += diff * diff;
    }
    return sum;
}

#[inline]
fn l1_dist(p1: &Point, p2: &Point) -> u64 {
    let mut sum: u64 = 0;
    for i in 0..DIM {
        sum += (p1[i] as i64 - p2[i] as i64).abs() as u64;
    }
    return sum;
}

#[inline]
fn get_position(p: &Point, source: &Point) -> usize {
    let mut pos: usize = 0;
    for i in 0..DIM {
        if p[i] > source[i] {
            pos += 1 << i;
        }
    }
    return pos;
}
#[inline]
fn intersection(p: &Point, metric: u32) -> Vec<Point> {
    if DIM != 2 {
        panic!("DIM should be 2");
    }
    let mut result: Vec<Point> = Vec::with_capacity(BLK_CELLS);
    let blk = block(p, SIDE_LEN, R);
    let mut cross_point: Point = [0u64; DIM];
    for j in 0..DIM {
        cross_point[j] = blk[j] * SIDE_LEN + SIDE_LEN;
    }
    let dist;
    if metric == 2 {
        dist = l2_dist(p, &cross_point);
    } else if metric == 1 {
        dist = l1_dist(p, &cross_point);
    } else {
        panic!("metric should be L1 or L2");
    }
    let pos_ind = get_position(&cross_point, p);
    for i in 0..BLK_CELLS {
        let mut tem: Point = [0u64; DIM];
        let r_lp = if metric == 2 { R_L2 } else { R };
        if (dist > r_lp) && (i == pos_ind) {
            continue;
        }
        for j in 0..DIM {
            if (i >> j) & 1 == 1 {
                tem[j] = blk[j] + 1;
            } else {
                tem[j] = blk[j];
            }
        }
        result.push(tem);
    }
    return result;
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
    pub fn new(num_item: u64, apart: bool) -> Self {
        // generate random sk and pk
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng);
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE;
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(DIM);
        let mut _okvsgen: Vec<okvs::OkvsGen> = Vec::with_capacity(DIM);
        let window: usize = if apart == true {
            BLK_CELLS * (2 * R + 1) as usize
        } else {
            (2 * R + 1) as usize
        };
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
                // for each possible value in [2R+1]
                let min = pt[i] - R as u64;
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

    pub fn msg_apart(&mut self, pt_set: &Vec<Point>) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(DIM);
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        let possible_vals = 2 * R as usize + 1;
        // for each dimension
        for i in 0..DIM {
            // for each receiver's point pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            {
                let blk = block(pt, SIDE_LEN, R);
                // for each possible value in [2R+1]*BLK_CELLS
                let mut cel: Point = [0u64; DIM];
                for k in 0..BLK_CELLS {
                    for j in 0..DIM {
                        if (k >> j) & 1 == 1 {
                            cel[j] = blk[j] + 1;
                        } else {
                            cel[j] = blk[j];
                        }
                    }
                    let key = hash64(&cel);
                    let min = pt[i] - R as u64;
                    for j in 0..possible_vals {
                        let key_ij = hash64(&(min + j as u64));
                        list.push((key ^ key_ij, pre_window[k * possible_vals + j]));
                    }
                }
            }
            result.push(self._okvsgen[i].encode(&list));
            list.clear();
        }
        self._okvsgen.clear();
        self._pre_data.clear();
        return result;
    }

    pub fn lp_msg_apart(&mut self, pt_set: &Vec<Point>, metric: u32) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(DIM);
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        let possible_vals = 2 * R as usize + 1;
        // for each dimension
        for i in 0..DIM {
            // for each receiver's point pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            {
                // for each possible value in [2R+1]*BLK_CELLS
                let cels = intersection(pt, metric);
                for k in 0..BLK_CELLS {
                    let key;
                    if k >= cels.len() {
                        key = rand::random::<u64>();
                    } else {
                        key = hash64(&cels[k]);
                    }
                    let min = pt[i] - R as u64;
                    for j in 0..possible_vals {
                        let key_ij = hash64(&(min + j as u64));
                        let tem = pre_window[k * possible_vals + j];
                        let mut diff_abs = if j as u64 > R {
                            j as u64 - R
                        } else {
                            R - j as u64
                        };
                        if metric == 2 {
                            diff_abs *= diff_abs;
                        }
                        list.push((key ^ key_ij, (tem.0, tem.1 + Scalar::from(diff_abs))));
                    }
                }
            }
            result.push(self._okvsgen[i].encode(&list));
            list.clear();
        }
        self._okvsgen.clear();
        self._pre_data.clear();
        return result;
    }

    #[inline]
    pub fn post_process(&mut self, msg_sender: &okvs::Encoding) -> u32 {
        for (u, v) in msg_sender.iter() {
            if (self.sk * u) == *v {
                return 1;
            }
        }
        return 0;
    }
    #[inline]
    pub fn post_process_apart(&mut self, msg_sender: &okvs::PointPair) -> u32 {
        if (self.sk * msg_sender.0) == msg_sender.1 {
            return 1;
        }
        return 0;
    }
    #[inline]
    pub fn lp_post_process_apart(&mut self, msg_sender: &(okvs::PointPair, HashSet<u32>)) -> u32 {
        let x = hash64(
            &(msg_sender.0 .1 - self.sk * msg_sender.0 .0)
                .compress()
                .to_bytes(),
        ) as u32;
        if msg_sender.1.contains(&x) {
            return 1;
        }
        return 0;
    }
    pub fn output(&mut self, msg_sender: &okvs::Encoding, window: usize) -> u32 {
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
    pub fn output_apart(&mut self, msg_sender: &okvs::Encoding) -> u32 {
        let mut count = 0;
        for (u, v) in msg_sender.iter() {
            if (self.sk * u) == *v {
                count += 1;
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
    _coins_lp_set: Vec<HashSet<u32>>,
    _coins_lp: Vec<RistrettoPoint>,
}
impl Sender {
    pub fn new(num_item: u64, pk_rec: RistrettoPoint, apart: bool, metric: u32) -> Self {
        let window: usize = if apart == true { 1 } else { BLK_CELLS };
        let m = num_item * window as u64;
        let mut rng = rand::thread_rng();
        let pk = pk_rec;
        let mut _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)> =
            Vec::with_capacity(m as usize);
        let metric_window = if metric == 2 {
            R_L2 as usize + 1
        } else {
            R as usize + 1
        };
        let mut _coins_lp_set: Vec<HashSet<u32>> = Vec::with_capacity(m as usize);
        let mut _coins_lp: Vec<RistrettoPoint> = Vec::with_capacity(m as usize);
        for _ in 0..m {
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);
            let c = Scalar::random(&mut rng);
            _coins.push((&a * RISTRETTO_BASEPOINT_TABLE, a * pk_rec, b));
            _coins_lp.push(&c * RISTRETTO_BASEPOINT_TABLE);
            let mut hashtab: HashSet<u32> = HashSet::with_capacity(metric_window);
            for i in 0..metric_window as u64 {
                let g_j = &(c + b * Scalar::from(i)) * RISTRETTO_BASEPOINT_TABLE;
                hashtab.insert(hash64(&g_j.compress().to_bytes()) as u32);
            }
            _coins_lp_set.push(hashtab);
        }
        return Sender {
            m,
            window,
            pk,
            _coins,
            _coins_lp_set,
            _coins_lp,
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

    #[inline]
    pub fn send_msg_single(
        &self,
        encodings: &Vec<okvs::Encoding>,
        pt: &Point,
        index: usize,
    ) -> okvs::Encoding {
        let mut blk: Point = [0u64; DIM];
        let coin_window = &self._coins[index..index + self.window];
        let mut uv: okvs::Encoding =
            vec![(RistrettoPoint::identity(), RistrettoPoint::identity()); BLK_CELLS];
        let mut tem: okvs::PointPair;
        let cel = cell(pt, SIDE_LEN);
        // for each possible block
        for (i, coins) in coin_window.iter().enumerate() {
            for j in 0..DIM {
                if (i >> j) & 1 == 1 {
                    blk[j] = cel[j] - 1;
                } else {
                    blk[j] = cel[j];
                }
            }
            let key = hash64(&blk);
            // for each dimension
            for j in 0..DIM {
                let key_ij = hash64(&(pt[j] as u64));
                tem = okvs::okvs_decode(&encodings[j], key ^ key_ij);
                uv[i].0 += tem.0;
                uv[i].1 += tem.1;
            }
            // finalize
            uv[i].0 = coins.2 * uv[i].0 + coins.0;
            uv[i].1 = coins.2 * uv[i].1 + coins.1;
        }
        return uv;
    }

    #[inline]
    pub fn send_msg_single_apart(
        &self,
        encodings: &Vec<okvs::Encoding>,
        pt: &Point,
        index: usize,
    ) -> okvs::PointPair {
        let coins = &self._coins[index];
        let mut uv: okvs::PointPair = (RistrettoPoint::identity(), RistrettoPoint::identity());
        let mut tem: okvs::PointPair;
        let cel = cell(pt, SIDE_LEN);
        let key = hash64(&cel);
        // for each dimension
        for j in 0..DIM {
            let key_ij = hash64(&(pt[j] as u64));
            tem = okvs::okvs_decode(&encodings[j], key ^ key_ij);
            uv.0 += tem.0;
            uv.1 += tem.1;
        }
        // finalize
        uv.0 = coins.2 * uv.0 + coins.0;
        uv.1 = coins.2 * uv.1 + coins.1;

        return uv;
    }

    #[inline]
    pub fn lp_send_msg_single_apart(
        &self,
        encodings: &Vec<okvs::Encoding>,
        pt: &Point,
        index: usize,
    ) -> (okvs::PointPair, HashSet<u32>) {
        let coins = &self._coins[index];
        let mut uv: okvs::PointPair = (RistrettoPoint::identity(), RistrettoPoint::identity());
        let mut tem: okvs::PointPair;
        let cel = cell(pt, SIDE_LEN);
        let key = hash64(&cel);
        // for each dimension
        for j in 0..DIM {
            let key_ij = hash64(&(pt[j] as u64));
            tem = okvs::okvs_decode(&encodings[j], key ^ key_ij);
            uv.0 += tem.0;
            uv.1 += tem.1;
        }
        // finalize
        uv.0 = coins.2 * uv.0 + coins.0;
        uv.1 = coins.2 * uv.1 + coins.1 + &self._coins_lp[index];

        return (uv, self._coins_lp_set[index].clone());
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

    pub fn msg_apart(
        &mut self,
        encodings: &Vec<okvs::Encoding>,
        pt_set: &Vec<Point>,
    ) -> okvs::Encoding {
        let mut uv: okvs::Encoding =
            vec![(RistrettoPoint::identity(), RistrettoPoint::identity()); self.m as usize];
        let mut tem: (RistrettoPoint, RistrettoPoint);
        // for each senders's point pt
        for (i, (pt, coins)) in pt_set.iter().zip(self._coins.iter()).enumerate() {
            let cel = cell(pt, SIDE_LEN);
            let key = hash64(&cel);
            // for each dimension
            for j in 0..DIM {
                let key_ij = hash64(&(pt[j] as u64));
                tem = okvs::okvs_decode(&encodings[j], key ^ key_ij);
                uv[i].0 += tem.0;
                uv[i].1 += tem.1;
            }
            // finalize
            uv[i].0 = coins.2 * uv[i].0 + coins.0;
            uv[i].1 = coins.2 * uv[i].1 + coins.1;
        }
        self._coins.clear();
        return uv;
    }
}
