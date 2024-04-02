use std::{
    fmt,
    ops::{Index, IndexMut}
};

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SquareMatrix<const S: usize> {
    Null,
    NonNull([[f64; S]; S])
}

impl<const S: usize> SquareMatrix<S> {
    #[inline]
    pub fn from(matrix: [[f64; S]; S]) -> Self {
        Self::NonNull(matrix)
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        match self {
            SquareMatrix::Null       => 0.,
            SquareMatrix::NonNull(m) => m[i][j]
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut r = Vec::<u8>::with_capacity(S*S*8+1);
        match self {
            Self::Null => {
                r.resize(S*S*8+1, 0);
            },
            Self::NonNull(m) => {
                for i in 0..S {
                    for j in 0..S {
                        for b in m[i][j].to_le_bytes() {
                            r.push(b);
                        }
                    }
                }
                r.push(1);
            }
        }
        r
    }
    pub fn from_bytes(matrix: &[u8]) -> Result<Self, String> {
        if matrix.len() != S*S*8+1 {
            return Err(format!("wrong size ({}), expected {}", matrix.len(), S*S*8+1))
        }

        match matrix.last() {
            Some(&0) => Ok(Self::Null),
            Some(&1) => {
                let mut m = [[0_f64; S];S];
                for (idx, number) in matrix.chunks_exact(8).enumerate() {
                    let (i, j) = (idx / S, idx % S);
                    m[i][j] = f64::from_le_bytes(number.try_into().unwrap());
                }
                Ok(Self::from(m))
            },
            Some(&e) => Err(format!("Enum variant not recognized: {}", e)),
            None => Err("Empty stream".to_string())
        }
    }
}

impl<const S: usize> Index<(usize, usize)> for SquareMatrix<S> {
    type Output = f64;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        match self {
            SquareMatrix::Null       => &0_f64,
            SquareMatrix::NonNull(m) => &m[i][j]
        }
    }
}

impl<const S: usize> fmt::Display for SquareMatrix<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SquareMatrix::Null       => write!(f, "Null matrix")?,
            SquareMatrix::NonNull(m) => {
                for i in 0..S {
                    for j in 0..S {
                        write!(f, "\t{}", m[i][j])?;
                    }
                    write!(f, "\n")?;
                }
            }
        }
        write!(f, "")
    }
}


#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vector<const S: usize> {
    NonNull([f64; S])
}

impl<const S: usize> Vector<S> {
    #[inline]
    pub fn new() -> Self {
        Vector::NonNull([0_f64; S])
    }

    #[inline]
    pub fn get(&self, i: usize) -> f64 {
        match self {
            Vector::NonNull(v) => v[i]
        }
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut r = Vec::<u8>::with_capacity(S*8+1);
        match self {
            Self::NonNull(v) => {
                for &i in v {
                    for b in i.to_le_bytes() {
                        r.push(b);
                    }
                }
                r.push(1);
            }
        }
        r
    }
    pub fn from_bytes(vector: &[u8]) -> Result<Self, String> {
        if vector.len() != S*8+1 {
            return Err(format!("wrong size ({}), expected {}", vector.len(), S*8+1))
        }

        match vector.last() {
            Some(&1) => {
                let mut v = [0_f64; S];
                for (i, number) in vector.chunks_exact(8).enumerate() {
                    v[i] = f64::from_le_bytes(number.try_into().unwrap());
                }
                Ok(Self::from(v))
            },
            Some(&e) => Err(format!("Enum variant not recognized: {}", e)),
            None => Err("Empty stream".to_string())
        }
    }
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut f64> {
        match self {
            Vector::NonNull(v) => v.iter_mut()
        }
    }
    pub fn iter(&self) -> impl Iterator<Item=&f64> {
        match self {
            Vector::NonNull(v) => v.iter()
        }
    }
    pub fn to_vec(&self) -> Vec<f64> {
        match self {
            Vector::NonNull(v) => v.to_vec()
        }
    }
    pub fn from_vec(vec: &Vec<f64>) -> Self {
        let mut v = Self::new();
        for (i, &value) in vec.iter().enumerate() {
            v[i] = value;
        }
        v
    }
}

impl<const S: usize> Index<usize> for Vector<S> {
    type Output = f64;
    #[inline]
    fn index(&self, i: usize) -> &f64 {
        match self {
            Self::NonNull(v) => &v[i]
        }
    }
}

impl<const S: usize> IndexMut<usize> for Vector<S> {
    fn index_mut(&mut self, i: usize) -> &mut f64 {
        match self {
            Self::NonNull(v) => &mut v[i]
        }
    }
}

impl<const S: usize> IntoIterator for Vector<S> {
    type Item = f64;
    type IntoIter = std::array::IntoIter<Self::Item, S>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::NonNull(v) => v.into_iter()
        }
    }
}


impl<const S: usize> fmt::Display for Vector<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Vector::NonNull(v) => {
                for i in 0..S {
                    write!(f, "\t{}", v[i])?;
                }
            }
        }
        write!(f, "")
    }
}


impl<const S: usize> ::std::convert::From<[f64; S]> for Vector<S> {
    #[inline]
    fn from(vector: [f64; S]) -> Self {
        Self::NonNull(vector)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn save_load() {
        let a = SquareMatrix::<3>::Null;
        assert_eq!(a, SquareMatrix::<3>::from_bytes(&a.to_bytes()).unwrap());

        let b = SquareMatrix::<3>::from([[1_f64; 3]; 3]);
        assert_eq!(b, SquareMatrix::<3>::from_bytes(&b.to_bytes()).unwrap());

        let b = Vector::<3>::from([1_f64; 3]);
        assert_eq!(b, Vector::<3>::from_bytes(&b.to_bytes()).unwrap());

        let b = Vector::<3>::from([1_f64; 3]);
        let mut b_mod = b.to_bytes();
        b_mod[7] = 1;
        assert_ne!(b, Vector::<3>::from_bytes(&b_mod).unwrap());
    }
}
