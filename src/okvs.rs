const EPSILON: f64 = 0.25;
const STAT_BITS: u64 = 40;
const FACTOR: f64 = STAT_BITS as f64 * 1.44 as f64;
const HASH_SEED: u64 = 0x1234567890abcdef;
const KAPPA: u64 = STAT_BITS;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use fxhash::hash64;
use rand::rngs::OsRng;

type Point = RistrettoPoint;

pub struct GBF {
    data: Vec<Point>,
    m: u64,
    n: u64,
}

impl GBF {
    pub fn new(num_item: u64) -> Self {
        let m: u64 = (num_item as f64 * FACTOR).floor() as u64;
        let data: Vec<Point> = (0..m).map(|_| Point::random(&mut OsRng)).collect();
        return GBF {
            data,
            m,
            n: num_item,
        };
    }
    pub fn num_items(&self) -> u64 {
        return self.n;
    }
    pub fn len(&self) -> u64 {
        return self.m;
    }
    pub fn num_hashes(&self) -> u64 {
        return KAPPA;
    }
    pub fn encode(&mut self, list: &Vec<(u64, Point)>) {
        let mut touched: Vec<bool> = vec![false; self.m as usize];
        let mut hashkey: [u64; STAT_BITS as usize] = [0u64; STAT_BITS as usize];
        for (key, value) in list.iter() {
            let seed = hash64(&(*key ^ HASH_SEED));
            for i in 0..STAT_BITS {
                hashkey[i as usize] = hash64(&(seed ^ i)) % self.m;
            }
            let mut temp = *value;
            let mut pos = 0usize;
            for i in 0..STAT_BITS as usize {
                let j = hashkey[i] as usize;
                if false == touched[j] {
                    pos = j;
                    touched[j] = true;
                }
                temp -= self.data[j];
            }
            assert!(pos != 0usize);
            self.data[pos] += temp;
        }
    }

    pub fn decode(&mut self, key: u64) -> Point {
        let seed: u64 = hash64(&(key ^ HASH_SEED));
        let mut hashkey = [0u64; STAT_BITS as usize];
        for i in 0..STAT_BITS {
            hashkey[i as usize] = hash64(&(seed ^ i)) % self.m;
        }

        let mut result = Point::identity();
        for i in 0..STAT_BITS as usize {
            result += self.data[hashkey[i] as usize];
        }
        return result;
    }
}

pub fn okvs_decode(data: &Vec<(Point, Point)>, key: u64) -> (Point, Point) {
    let pos_band_range: u64 = data.len() as u64 - KAPPA;
    let seed: u64 = hash64(&(key ^ HASH_SEED));
    let pos: usize = (seed % pos_band_range) as usize;
    let band: Vec<u64> = (0..KAPPA).map(|i| (seed >> i) & 0x01).collect();

    let mut result = (Point::identity(), Point::identity());

    for i in pos..KAPPA as usize + pos {
        if 0 == band[i - pos] {
            continue;
        }
        result.0 += data[i].0;
        result.1 += data[i].1;
    }
    return result;
}
pub struct OKVS {
    // data: Vec<(Point, Point)>,
    _data: Vec<(Scalar, Scalar)>,
    _matrix: Vec<(usize, usize, (Vec<Scalar>, (Scalar, Scalar)))>,
    m: u64,
    n: u64,
    epsilon: f64,
}
impl OKVS {
    pub fn new(num_item: u64) -> Self {
        let m: u64 = (num_item as f64 * (1.0 + EPSILON)).ceil() as u64;
        let _data: Vec<(Scalar, Scalar)> = (0..m)
            // .map(|_| (Scalar::ZERO, Scalar::ZERO))
            .map(|_| (Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)))
            .collect();
        return OKVS {
            // data: Vec::with_capacity(m as usize),
            _data,
            _matrix: Vec::with_capacity(num_item as usize),
            m,
            n: num_item,
            epsilon: EPSILON,
        };
    }
    pub fn num_items(&self) -> u64 {
        return self.n;
    }
    pub fn len(&self) -> u64 {
        return self.m;
    }
    pub fn expansion_rate(&self) -> f64 {
        return self.epsilon;
    }

    pub fn refresh(&mut self) {
        self._matrix.clear();
        self._data.clear();
        for _ in 0..self.m {
            self._data
                .push((Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)));
        }
    }

    pub fn encode(&mut self, list: &Vec<(u64, (Scalar, Scalar))>) -> Vec<(Point, Point)> {
        let mut data: Vec<(Point, Point)> = Vec::with_capacity(self.m as usize);
        let pos_band_range: u64 = self.m - KAPPA;
        for (key, value) in list.iter() {
            let seed = hash64(&(*key ^ HASH_SEED));
            let hash_pos = (seed % pos_band_range) as usize;
            let hash_band: Vec<Scalar> = (0..KAPPA)
                .map(|i| Scalar::from((seed >> i) & 0x01))
                .collect();

            self._matrix.push((hash_pos, 0usize, (hash_band, *value)));
        }
        // sort by position
        self._matrix.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        let mut pivots: usize = 0;
        for row in 0..self.n as usize {
            let (top, bot) = self._matrix.split_at_mut(row + 1);
            if let Some((pos, piv, (band, value))) = top.last_mut() {
                for i in 0..KAPPA as usize {
                    if Scalar::ZERO == band[i] {
                        continue;
                    }
                    *piv = *pos + i;
                    pivots += 1;
                    // if the pivot is not one, divide the row to normalize it
                    if Scalar::ONE != band[i] {
                        let inv = band[i].invert();
                        band[i] = Scalar::ONE;
                        for j in i + 1..KAPPA as usize {
                            band[j] *= inv;
                        }
                        value.0 *= inv;
                        value.1 *= inv;
                    }
                    // update the band from the following rows
                    // for j in row + 1..self.n as usize {
                    for (pos_j, piv_j, (band_j, value_j)) in bot.iter_mut() {
                        // skip if the position is out of range
                        if *pos_j > *piv {
                            continue;
                        }
                        // if found the non-zero position, subtract the band
                        *piv_j = *piv - *pos_j;
                        let multplier = band_j[*piv_j];
                        if Scalar::ZERO != multplier {
                            if Scalar::ONE != multplier {
                                for k in i + 1..KAPPA as usize {
                                    band_j[*piv_j + k - i] -= band[k] * multplier;
                                }
                                value_j.0 -= value.0 * multplier;
                                value_j.1 -= value.1 * multplier;
                            } else {
                                for k in i + 1..KAPPA as usize {
                                    band_j[*piv_j + k - i] -= band[k];
                                }
                                value_j.0 -= value.0;
                                value_j.1 -= value.1;
                            }
                        }
                        band_j[*piv_j] = Scalar::ZERO;
                    }
                    break;
                }
            }
        }
        // band should not be all-zero
        assert_eq!(pivots, self.n as usize);
        // back substitution
        for (pos, piv, (band, value)) in self._matrix.iter().rev() {
            let (mut val, mut val2) = *value;
            for i in 0..KAPPA as usize {
                if Scalar::ZERO == band[i] {
                    continue;
                }
                if (i + *pos) == *piv {
                    continue;
                }
                val -= self._data[i + *pos].0 * band[i];
                val2 -= self._data[i + *pos].1 * band[i];
            }
            self._data[*piv] = (val, val2);
        }
        // finalize
        for i in 0..self.m as usize {
            data.push((
                &self._data[i].0 * RISTRETTO_BASEPOINT_TABLE,
                &self._data[i].1 * RISTRETTO_BASEPOINT_TABLE,
            ));
        }
        self._matrix.clear();
        // self._data.clear();
        return data;
    }
}

#[test]
fn okvs_test() {
    use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT;
    use std::time::Instant;

    let n: u64 = 2000;
    let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
    for j in 0..n {
        list.push((j + 1, (Scalar::ONE, Scalar::ONE)));
    }
    let mut okvsmod = OKVS::new(n);

    let now = Instant::now();

    let data = okvsmod.encode(&list);

    let elapsed = now.elapsed();
    println!(
        "{} items, Elapsed Time for Encoding (optimize=0): {:.2?}",
        n, elapsed
    );
    for i in 0..n {
        let decoding = okvs_decode(&data, i + 1);
        assert_eq!(
            decoding,
            (RISTRETTO_BASEPOINT_POINT, RISTRETTO_BASEPOINT_POINT)
        );
    }
}
