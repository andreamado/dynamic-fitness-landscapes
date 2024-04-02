use std::{
    ops::Index,
    fmt::{self, Write}
};
use rand::Rng;

/// Returns the number of genotypes in a landscape with L biallelic loci
pub const fn landscape_size<const L: usize>() -> usize {
    2_usize.pow(L as u32)
}

/// Returns an vector with all the possible binary sequences of length L
pub fn possible_sequences<const L: usize>() -> Vec<[u8; L]> {
    let mut seq_list = Vec::<[u8; L]>::with_capacity(landscape_size::<L>());
    for g in 0..landscape_size::<L>() {
        let mut seq = [0u8; L];
        let mut t = g;
        for i in 0..L {
            seq[i] = (t % 2) as u8;
            t /= 2;
        }
        seq_list.push(seq);
    }
    seq_list
}

/// Type that represents a genotype.
/// It is meant to abstract out genotype representation so easy alternative implementations can be
/// supplied. The sequence is supposed to be accessed only through the provided methods.
#[derive(Copy, Clone, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
pub struct Genotype<const L: usize> {
    seq: [u8; L]
}

impl<const L: usize> Genotype<L> {
    /// Creates a new genotype with all alleles set to zero.
    #[inline]
    pub const fn new() -> Self {
        Genotype {
            seq: [0; L]
        }
    }

    /// Creates a new genotype from the given sequence.
    #[inline]
    pub fn from_sequence(sequence: &[u8]) -> Self {
        Genotype {
            seq: sequence.try_into().expect("slice lenght is wrong")
        }
    }

    /// Creates a new genotype with all alleles set to zero.
    #[inline]
    pub fn random() -> Self {
        let mut rng = rand::thread_rng();
        let side = rand::distributions::Uniform::new(0u8, 2u8);

        let mut seq = [0u8; L];
        for i in seq.iter_mut() {
            *i = rng.sample(side);
        }
        Genotype { seq }
    }

    pub fn from_index(index: usize) -> Self {
        let mut seq = [0u8; L];
        for i in 0..L {
            seq[i] = ((index / 2_usize.pow(i as u32)) % 2) as u8;
        }
        Genotype { seq }
    }

    /// Switches the ith allele of the genotype
    #[inline]
    pub fn mutate(&mut self, i: usize) -> &Self {
        self.seq[i] = 1 - self.seq[i];
        self
    }

    /// Clones the genotype and mutates the ith allele
    #[inline]
    pub fn cmutate(self, i: usize) -> Self {
        *self.clone().mutate(i)
    }

    #[inline]
    pub fn iter(&self) -> std::slice::Iter<u8> {
        self.seq.iter()
    }

    #[inline]
    pub fn to_vec(&self) -> Vec<u8> {
        self.seq.iter().map(|&x| x).collect()
    }

    pub fn sum(&self) -> usize {
        self.seq.iter().sum::<u8>() as usize
    }

    pub fn n_differences(&self, g2: &Self) -> usize {
        self.iter().zip(g2.iter()).map(|(l1, l2)| (*l1 as i16 - *l2 as i16).abs() as usize).sum()
    }

    pub fn index(&self) -> usize {
        self.iter().enumerate().fold(0, |acc, (i, s)| acc + 2_usize.pow(i as u32)*(*s as usize))
    }
    pub fn order(&self) -> usize {
        self.index() + 2_usize.pow(L as u32) * self.sum()
    }
}

impl<const L: usize> Index<usize> for Genotype<L> {
    type Output = u8;
    #[inline]
    fn index(&self, i: usize) -> &Self::Output {
        &self.seq[i]
    }
}

impl<const L: usize> IntoIterator for Genotype<L> {
    type Item = u8;
    type IntoIter = std::array::IntoIter<Self::Item, L>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.seq.into_iter()
    }
}

impl<const L: usize> fmt::Display for Genotype<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut display = String::new();
        for allele in self.seq {
            write!(display, " {}", allele)?
        }
        write!(f, "{}", &display[1..])
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general() {
        let genotype1 = Genotype::<5>::from_sequence(&[0, 0, 0, 1, 0]);
        let mut genotype2 = Genotype::<5>::new();
        genotype2.mutate(3);

        assert_eq!(genotype1, genotype2);
    }
}
