import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

def con_shift(array,x,n =2):
    shift = int(len(x) / n)
    N = len(x)
    kvec = np.arange(len(x))
    J = np.complex(0, 1)
    g = signal.unit_impulse(N,shift-1)
    g_fft = np.fft.fft(g)
    yft_shift = np.fft.fft(y)
    # de = np.exp(-2 * np.pi * J * kvec * shift / len(x))

    y_shift = np.real(np.fft.ifft(yft_shift*g_fft))
    return y_shift
x = np.arange(-10,10,.05)
y = np.exp(-0.5*x**2/(0.8**2))

gauss_shift = con_shift(y,x,3)


# plt.plot(x,y)
# plt.plot(x,gauss_shift)
# plt.title('1/3 shift of Gaussian')
# plt.savefig("gauss_shift_new.png",dpi = 300, bbox_inches = 'tight')
# plt.show()

def correl(a1, a2):
    a1_ft = np.fft.fft(a1)
    a2_ft = np.fft.fft(a2)
    cor = np.fft.ifft(a1_ft*np.conj(a2_ft))
    return cor

cor = correl(y,y)
# plt.plot(x,cor)
# plt.title('Correlation of two Gaussian')
# plt.savefig("Corr_Gauss.png",dpi = 300, bbox_inches = 'tight')
# plt.show()

def shift_cor(a1,a2, n):
    new_a2 = con_shift(a2, x, n)
    new_a2 = new_a2/new_a2.sum()
    a1 = a1/a1.sum()
    output = correl(a1,new_a2)
    return output,new_a2,a1
n = 3
Shift_Gauss, shift_fun , new_y= shift_cor(y,y,n)
# plt.plot(x,shift_fun,'r',label = 'Shifted Gauss')
# plt.plot(x,new_y,'g',label = "Original Gauss")
# plt.plot(x,Shift_Gauss,label = 'Correlated Function')
# plt.title(f'Correlation of Gaussian with shift of {1/n} of its length')
# plt.legend()
# plt.savefig(f'Corr_Shift_{1/n}.png',dpi = 300, bbox_inches = 'tight')
# plt.show()


def conv_safe(f,g):
    N = len(f) + len(g) - 1
    N1 = N - len(f)
    N2 = N - len(g)
    f_new = np.pad(f, (0, N1))
    g_new = np.pad(g, (0, N2))
    fg = np.abs(np.fft.ifft(np.fft.fft(f_new)*np.fft.fft(g_new)))
    return fg
xg = np.arange(-10,10,.01)
g = np.sin(xg)
ygconv = conv_safe(y,g)
plt.plot(ygconv,'orange')
plt.title('Convolution of Gauss and Sin Functions')
plt.savefig("Conv_safe_Gauss_Sin.png",dpi = 300, bbox_inches = 'tight')
plt.show()