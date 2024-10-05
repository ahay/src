k = var('k')
def ost(t, w):
    return plot(w*cos(t) - k*sin(t), k, -2, 2)
def env(w):
    return plot(sqrt(w^2 + k^2), k, -2, 2, thickness=5)
lines = [ost(t*pi/20,1) for t in range(-10, 11)]
plot = sum(lines)+env(1) 
plot.save(frame=True, axes=False, axes_labels=['Wavenumber (rad)', 'Frequency (rad)'], filename='junk_sage.pdf')
