import numpy as np

t = np.arange(0.5, 3.2, 0.001)

# smooth ramp (example: cosine shape)
a = 0.75 * 0.5 * (1 - np.cos(np.pi * (t - 0.5) / (2.0 - 0.5)))


print("(0.00   0.00) \n")
print("(0.20   0.00) \n")
print("(0.30   0.00) \n")
print("(0.40   0.00) \n")

for ti, ai in zip(t, a):
    print(f"({ti:.4f}   {ai:.6f})") 


print("(3.80   0.00) \n")
print("(4.20   0.00) \n")
print("(4.70   0.00) \n")
print("(5.00   0.00) \n")

