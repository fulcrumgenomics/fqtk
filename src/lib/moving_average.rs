/// A simple moving average calculator.
/// Only requires that T is convertable to f64.
/// Uses space of window * size_of(T) bytes.
pub(crate) struct MovingAverage<T> {
    window: usize,
    values: Vec<T>,
    sum: f64,
    idx: usize,
    count: usize,
}

impl<T: Copy + Default + std::convert::Into<f64>> MovingAverage<T> {
    /// create a new moving average calculator with a window of `window` values.
    pub fn new(window: usize) -> Self {
        Self { window, values: vec![T::default(); window], sum: 0.0, idx: 0, count: 0 }
    }

    /// push a new value into the moving average calculator and get the new mean.
    pub fn push(&mut self, value: T) -> f64 {
        let old_value = self.values[self.idx];
        self.values[self.idx] = value;
        self.sum = self.sum + value.into() - old_value.into();
        self.idx = (self.idx + 1) % self.window;
        self.count += 1;
        self.mean()
    }

    /// get the current mean.
    #[inline]
    pub fn mean(&self) -> f64 {
        self.sum / (self.count.min(self.window) as f64)
    }
}

// write some tests for the calculator
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_moving_average() {
        let window_size = 3;
        let mut ma = MovingAverage::new(window_size);
        // NOTE the first value is always the mean
        // we use min of values added and window size to calculate the mean
        assert_eq!(ma.push(1), 1 as f64 / 1 as f64);
        assert_eq!(ma.push(2), (1 + 2) as f64 / 2 as f64);
        assert_eq!(ma.push(3), (1 + 2 + 3) as f64 / window_size as f64);
        assert_eq!(ma.push(4), (2 + 3 + 4) as f64 / window_size as f64);
        assert_eq!(ma.push(5), (3 + 4 + 5) as f64 / window_size as f64);
        assert_eq!(ma.push(6), (4 + 5 + 6) as f64 / window_size as f64);
    }
}
