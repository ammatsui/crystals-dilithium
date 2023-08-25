use crate::{params::*, sample::*, poly::*, packing::*};
use rand::Rng;

/* CPAPKE procedures */

pub fn keyGen() -> ([u8; PK_BYTES], [u8; SK_BYTES])
{
    let mut sk = [0u8; SK_BYTES];
    let mut pk = [0u8; PK_BYTES];

    let d  = rand::thread_rng().gen::<[u8; 32]>();
    
    let (mut rho, mut sigma) = ([0u8; 32], [0u8; 32]);

    let mut tmp = [0u8; 64];
    G(&d, &mut tmp);
    rho.copy_from_slice(&tmp[..32]);
    sigma.copy_from_slice(&tmp[32..]);

    let A = get_matrix(&rho);
    let mut s = get_noise(&sigma);
    let mut e = get_error(&sigma);

    s = Ntt(&s);
    e = Ntt(&e);

    let mut t: VecPoly<K> = m_mult_v(&A, &Ntt(&s));
    t = v_add(&t, &e);

    pack_pk(&mut pk, &t, &rho);
    pack_sk(&mut sk, &s);

    (pk, sk)
}