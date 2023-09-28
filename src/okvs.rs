
pub const EPSILON: f64 = 0.15;
pub const STAT_BITS: u64 = 40;
const STAT_MASK: u64 = (1 << STAT_BITS) - 1;
const POS_MASK: u64 = 1 << STAT_BITS;

use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use fxhash::hash64;
use rand::rngs::OsRng;
use std::collections::HashMap;

pub fn okvs_encode(list: &HashMap<u64, RistrettoPoint>) -> Vec<RistrettoPoint> {
    // Implementation goes here
    let n: u64 = list.len() as u64;
    let m: u64 = (n as f64 * (1.0 + EPSILON)).ceil() as u64;
    let mut pivots: Vec<u64> = vec![0u64; n as usize];
    let pos_band_range: u64 = m - STAT_BITS;
    let mut matrix: Vec<(u64, (u64, RistrettoPoint))> = Vec::with_capacity(n as usize);
    let mut result: Vec<RistrettoPoint> =
        (0..m).map(|_| RistrettoPoint::random(&mut OsRng)).collect();
    let mut hash_pos: u64;
    let mut hash_band: u64;

    for (&key, &value) in list.iter() {
        hash_pos = hash64(&key) % pos_band_range;
        hash_band = hash64(&hash_pos) & STAT_MASK;
        matrix.push((hash_pos, (hash_band, value)));
    }
    // sort by position
    matrix.sort_by(|a, b| a.0.cmp(&b.0));
    for row in 0..n {
        let (pos, (band, value)) = matrix[row as usize];
        for i in 0..STAT_BITS {
            if 0 == (band & (POS_MASK >> i)) {
                continue;
            }
            pivots[row as usize] = i + pos;
            for j in row + 1..n {
                if matrix[j as usize].0 > pivots[row as usize] {
                    continue;
                }
                // shift bands to get bits under the pivot
                let to_compare: u64 =
                    matrix[j as usize].1 .0 << (pivots[row as usize] - matrix[j as usize].0);
                if 0 != (to_compare & POS_MASK) {
                    matrix[j as usize].1 .0 ^= band << (matrix[j as usize].0 - pos);
                    matrix[j as usize].1 .1 += value;
                }
            }
            break;
        }
    }
    // back substitution
    for row in (0..n).rev() {
        let (pos, (band, value)) = matrix[row as usize];
        result[pos as usize] = value;
        for i in 0..STAT_BITS {
            if 0 == (band & (POS_MASK >> i)) {
                continue;
            }
            if i + pos == pivots[row as usize] {
                continue;
            }
            result[pos as usize] = result[pos as usize] - result[(i + pos) as usize];
        }
    }
    return result;
}

pub fn okvs_decode(okvs: &Vec<RistrettoPoint>, key: u64) -> RistrettoPoint {
    // Implementation goes here
    let m: u64 = okvs.len() as u64;
    let pos_band_range: u64 = m - STAT_BITS;
    let pos: u64;
    let band: u64;
    pos = hash64(&key) % pos_band_range;
    band = hash64(&pos) & STAT_MASK;
    let mut result = RistrettoPoint::identity();

    for i in 0..STAT_BITS {
        if 0 == (band & (POS_MASK >> i)) {
            continue;
        }
        result += okvs[(i + pos) as usize];
    }
    return result;
}
