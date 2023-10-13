print("Tugas 3")
print ("Nadya Putri Arisni")
print("NRP = 5009211079")

import cmath

# FFT 1 dimensi
def fft(x):
    N = len(x)
    # Kasus dasar: Jika panjang input kurang dari atau sama dengan 1, kembalikan input tersebut
    if N <= 1:
        return x

    # Bagi input menjadi bagian genap dan ganjil
    genap = fft(x[0::2])
    ganjil = fft(x[1::2])

    # Hitung faktor Twiddle dan gabungkan hasilnya
    T = [cmath.exp(-2j * cmath.pi * k / N) * ganjil[k] for k in range(N // 2)]
    return [genap[k] + T[k] for k in range(N // 2)] + [genap[k] - T[k] for k in range(N // 2)]

# FFT 2D
def fft2d(matriks):
    baris = len(matriks)
    kolom = len(matriks[0])

   
    for i in range(baris):
        matriks[i] = fft(matriks[i])

    # Transposisi matriks untuk persiapan FFT kolom
    matriks_transpos = [[matriks[i][j] for i in range(baris)] for j in range(kolom)]


    for j in range(kolom):
        matriks_transpos[j] = fft(matriks_transpos[j])

    # Transposisi hasil kembali ke orientasi aslinya
    hasil_matriks = [[matriks_transpos[j][i] for j in range(kolom)] for i in range(baris)]

    return hasil_matriks

# Fungsi untuk menghitung MFCC dari matriks hasil FFT
def compute_mfcc(fft_result, num_mfcc_coeffs):
    # Mengambil log dari magnitude spektrum
    log_magnitude = [[abs(val) for val in row] for row in fft_result]
    
    def dct(x, type=2):
        N = len(x)
        y = [0] * N
        for k in range(N):
            y[k] = sum([x[n] * cmath.exp(-1j * cmath.pi * k * (2 * n + 1) / (2 * N)) for n in range(N)])
            if type == 2:
                y[k] *= 2
        return y
    
    mfcc = [dct(row)[:num_mfcc_coeffs] for row in log_magnitude]
    
    return mfcc

# Contoh penggunaan dengan sinyal step berbeda amplitudo
A = 1

# Sinyal dengan amplitudo 1/2A
matriks_input = [
    [1/(2*A), 1/(2*A), 1/(2*A), 1/(2*A)],
    [1/(2*A), 1/(2*A), 1/(2*A), 1/(2*A)],
    [1/(2*A), 1/(2*A), 1/(2*A), 1/(2*A)],
    [1/(2*A), 1/(2*A), 1/(2*A), 1/(2*A)]
]

# Hitung FFT 2D
hasil_fft2d = fft2d(matriks_input)

# Hitung MFCC dari hasil FFT 2D
num_mfcc_coeffs = 13  # Jumlah koefisien MFCC yang diinginkan
mfcc = compute_mfcc(hasil_fft2d, num_mfcc_coeffs)

# Tampilkan hasil MFCC
print("MFCC:")
for row in mfcc:
    print(row)
    
import matplotlib.pyplot as plt
import numpy as np

# Fungsi untuk plot spektrum MFCC
def plot_mfcc(mfcc):
    # Konversi data kompleks ke magnitudo
    mfcc_mag = np.abs(mfcc)

    plt.figure(figsize=(10, 5))
    plt.imshow(np.transpose(mfcc_mag), cmap='viridis', origin='lower', aspect='auto')
    plt.colorbar(format="%+2.0f dB")
    plt.title('MFCC')
    plt.xlabel('Frame')
    plt.ylabel('MFCC Coefficients')
    plt.show()

# Plot MFCC dari hasil perhitungan sebelumnya
plot_mfcc(mfcc)
