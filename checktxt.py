import numpy as np
import matplotlib.pyplot as plt

file_path = "build/"
# |%%--%%| <1hQv8kKawq|YFcZQ5BtEY>

Bz = []
for i in range(0, 2001, 20):
    file = file_path + f"Bz_{i}.txt"
    Bz.append(np.loadtxt(file))
print(np.shape(Bz))
# |%%--%%| <YFcZQ5BtEY|ikW2GcdhZD>

for i in range(0, 101, 10):
    plt.imshow(Bz[i])
    plt.colorbar()
    plt.show()

# |%%--%%| <ikW2GcdhZD|dAI0a1tWaF>

plt.plot(Bz[33][150][:])

# |%%--%%| <dAI0a1tWaF|NuGrbXSxiU>

Ex = []
for i in range(0, 1001, 10):
    file = file_path + f"Ex_{i}.txt"
    Ex.append(np.loadtxt(file))

for i in range(0, 102, 10):
    plt.imshow(Ex[i], vmax=1, vmin=-1)
    plt.colorbar()
    plt.show()

# |%%--%%| <NuGrbXSxiU|gVDk4Ghvix>


Ey = []
for i in range(0, 1001, 10):
    file = file_path + f"Ey_{i}.txt"
    Ey.append(np.loadtxt(file))

for i in range(0, 102, 10):
    plt.imshow(Ey[i], vmax=1, vmin=-1)
    plt.colorbar()
    plt.show()

