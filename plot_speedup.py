import matplotlib.pyplot as plt

# Given data points
data_points = [
    (5, 0.0044225),
    (50, 0.009344),
    (500, 0.58331825),
    (5000, 119.5131235),
    (50000, 1447.1241211667)
]

# Extract x and y values from data points
x_values = [point[0] for point in data_points]
y_values = [point[1] for point in data_points]

# Plot the points
plt.plot(x_values, y_values, marker='o', linestyle='-')

# Labeling axes
plt.xlabel('Array Size')
plt.ylabel('Speedup')

# Title of the plot
plt.title('Array Size vs Speed-Up')

# Show the plot
plt.show()
