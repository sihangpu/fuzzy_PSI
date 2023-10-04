use crate::psi;

use std::sync::mpsc as std_mpsc;
use std::thread;
// use tokio::sync::mpsc;

pub fn setup(num_n: usize, num_m: usize, apart: bool, metric: u32) -> (psi::Receiver, psi::Sender) {
    let psi_rec = psi::Receiver::new(num_n as u64, apart);
    let psi_sed = psi::Sender::new(num_m as u64, psi_rec.publish_pk(), apart, metric);
    return (psi_rec, psi_sed);
}

// #[tokio::main]
// pub async fn run_async(
//     mut psi_rec: psi::Receiver,
//     psi_sed: psi::Sender,
//     data_r: Vec<psi::Point>,
//     data_s: Vec<psi::Point>,
// ) {
//     // let (sender, mut receiver) = mpsc::channel::<[okvs::PointPair; psi::BLK_CELLS]>(data_s.len());
//     let (done_tx, mut done_rx) = mpsc::channel::<()>(1);
//     let (sender, mut receiver) = mpsc::unbounded_channel();

//     let msg1 = psi_rec.msg(&data_r);

//     tokio::spawn(async move {
//         for i in 0..data_s.len() {
//             let msg = psi_sed.send_msg_single(&msg1, &data_s[i], i * psi::BLK_CELLS);
//             if let Err(e) = sender.send(msg) {
//                 println!("Failed to send message: {:?}", e);
//             }
//         }
//     });

//     let mut count = 0u32;
//     // Receiver task.
//     tokio::spawn(async move {
//         while let Some(msg2) = receiver.recv().await {
//             count += psi_rec.post_process(&msg2);
//         }
//         println!("count: {}", count);
//         done_tx.send(()).await.expect("Failed to send done signal");
//     });

//     // Wait for the tasks to finish.
//     done_rx.recv().await.expect("Failed to receive done signal");
// }

pub fn run_standard(
    mut psi_rec: psi::Receiver,
    psi_sed: psi::Sender,
    data_r: Vec<psi::Point>,
    data_s: Vec<psi::Point>,
) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>();
    let (sender, receiver) = std_mpsc::channel();

    let msg1 = psi_rec.msg(&data_r);

    // Sender thread
    thread::spawn(move || {
        for i in 0..data_s.len() {
            let msg = psi_sed.send_msg_single(&msg1, &data_s[i], i * psi::BLK_CELLS);
            if sender.send(msg).is_err() {
                println!("Receiver has been dropped!");
                break;
            }
        }
    });

    // Receiver thread
    thread::spawn(move || {
        let mut count = 0u32;
        for msg2 in receiver.iter() {
            count += psi_rec.post_process(&msg2);
        }
        println!("count: {}", count);
        done_tx.send(()).expect("Failed to send done signal");
    });

    // Wait for the tasks to finish.
    done_rx.recv().expect("Failed to receive done signal");
}

pub fn run_standard_apart(
    mut psi_rec: psi::Receiver,
    psi_sed: psi::Sender,
    data_r: Vec<psi::Point>,
    data_s: Vec<psi::Point>,
) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>();
    let (sender, receiver) = std_mpsc::channel();

    let msg1 = psi_rec.msg_apart(&data_r);

    // Sender thread
    thread::spawn(move || {
        for i in 0..data_s.len() {
            let msg = psi_sed.send_msg_single_apart(&msg1, &data_s[i], i);
            if sender.send(msg).is_err() {
                println!("Receiver has been dropped!");
                break;
            }
        }
    });

    // Receiver thread
    thread::spawn(move || {
        let mut count = 0u32;
        for msg2 in receiver.iter() {
            count += psi_rec.post_process_apart(&msg2);
        }
        println!("count: {}", count);
        done_tx.send(()).expect("Failed to send done signal");
    });

    // Wait for the tasks to finish.
    done_rx.recv().expect("Failed to receive done signal");
}

pub fn run_standard_lp(
    mut psi_rec: psi::Receiver,
    psi_sed: psi::Sender,
    data_r: Vec<psi::Point>,
    data_s: Vec<psi::Point>,
    metric: u32,
) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>();
    let (sender, receiver) = std_mpsc::channel();

    let msg1 = psi_rec.lp_msg_apart(&data_r, metric);

    // Sender thread
    thread::spawn(move || {
        for i in 0..data_s.len() {
            let msg = psi_sed.lp_send_msg_single_apart(&msg1, &data_s[i], i);
            if sender.send(msg).is_err() {
                println!("Receiver has been dropped!");
                break;
            }
        }
    });

    // Receiver thread
    thread::spawn(move || {
        let mut count = 0u32;
        for msg2 in receiver.iter() {
            count += psi_rec.lp_post_process(&msg2);
        }
        println!("count: {}", count);
        done_tx.send(()).expect("Failed to send done signal");
    });

    // Wait for the tasks to finish.
    done_rx.recv().expect("Failed to receive done signal");
}
