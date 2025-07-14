import numpy as np
a = np.array([[296.62,0.011057], [299.62, 0.01443823], [293.62,0.007848]])
x= a[:, 0]
y= a[:, 1]
print(x)
print(y)
m, n,_ = np.polyfit(x, y, 2)
print(f"m={m}, n={n},min={-n/(2*m)}")    