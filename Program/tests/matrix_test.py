import numpy as np
import matplotlib.pyplot as plt

# Your list of values
values = np.array([10, 20, 30, 40, 50])

# Compute ratio matrix: row value ÷ column value
denom = np.where(values == 0, 1e-9, values)
ratio_matrix = values[:, None] / denom[None, :]

# Labels
value_labels = [str(v) for v in values]

# Plot: wide and short to squash vertically
fig, ax = plt.subplots(figsize=(7, 3))  # Width = 7, Height = 3 (squashed)

# Hide bottom x-axis, keep y-axis
ax.set_xticks([])
ax.set_yticks(np.arange(len(values)))
ax.set_yticklabels(value_labels)

# Top x-axis with column labels
ax_top = ax.secondary_xaxis('top')
ax_top.set_xticks(np.arange(len(values)))
ax_top.set_xticklabels(value_labels)
ax_top.set_xlabel('Compared *to* value (denominator)', fontsize=11)

# Annotate with 6 decimal places
for i in range(len(values)):
    for j in range(len(values)):
        ratio = ratio_matrix[i, j]
        ax.text(j, i, f'{ratio:.6f}', ha='center', va='center', fontsize=9)

# Grid lines
for i in range(len(values) + 1):
    ax.axhline(i - 0.5, color='black', linewidth=0.5)
    ax.axvline(i - 0.5, color='black', linewidth=0.5)

# Final touches
ax.set_ylabel('Original value (numerator)', fontsize=11)
ax.set_xlim(-0.5, len(values) - 0.5)
ax.set_ylim(len(values) - 0.5, -0.5)
ax.set_title('Ratio Matrix: (Row ÷ Column)', fontsize=13)
plt.tight_layout()
plt.show()

