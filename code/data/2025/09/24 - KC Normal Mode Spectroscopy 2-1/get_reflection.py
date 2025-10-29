import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data_from_fit.csv")

x = df["x"]
y = df["y"]

x_det_0 = x[34]
y_det_0 = y[34]

y_max_det = max(y)

print(f"Reflectivity at 0 detuning is {y_det_0 / y_max_det:.6f}")

plt.scatter(x[34], y[34], c="r")
plt.plot(x, y)
plt.title("Splitting @ 2-1'")
plt.show()
