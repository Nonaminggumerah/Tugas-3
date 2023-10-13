import cmath

# Handmade implementation of the Cooley-Tukey algorithm
# Make sure that the input signal is a power of two!

def omega(p, q):
    return cmath.exp((2.0 * cmath.pi * 1j * q) / p)

def fft1(signal):
    n = len(signal)
    if n == 1:
        return signal
    else:
        Feven = fft1([signal[i] for i in range(0, n, 2)])
        Fodd = fft1([signal[i] for i in range(1, n, 2)])
        combined = [0] * n
        for m in range(n // 2):
            combined[m] = Feven[m] + omega(n, -m) * Fodd[m]
            combined[m + n // 2] = Feven[m] - omega(n, -m) * Fodd[m]
        return combined

# Simple samples
s = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

fft_a = fft1(s)
fft_b = [float(x.real) for x in fft_a]

# Same
is_close = all(abs(a - b) < 1e-6 for a, b in zip(fft_a, fft_b))
assert is_close

# Two sine waves, 40 Hz + 90 Hz
t = [i / 500 for i in range(500)]
s = [cmath.sin(40 * 2 * cmath.pi * i) + 0.5 * cmath.sin(90 * 2 * cmath.pi * i) for i in t]

fft_a = fft1(s)
fft_b = [float(x.real) for x in fft_a]

# Check results
is_close = all(abs(a - b) < 1e-6 for a, b in zip(fft_a, fft_b))
assert is_close

print("Results match!")

