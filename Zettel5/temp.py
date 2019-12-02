import numpy as np
import matplotlib.pyplot as plt

u = np.array([[4], [2], [1]])


def step(x, rule_b):
    # The columns contains the L, C, R values
    # of all cells.
    y = np.vstack((np.roll(x, 1), x,
                   np.roll(x, -1))).astype(np.int8)
    # We get the LCR pattern numbers between 0 and 7.
    z = np.sum(y * u, axis=0).astype(np.int8)
    print(rule_b)
    # We get the patterns given by the rule.
    return rule_b[7 - z]


def generate(rule, size=7, steps=2):
    # Compute the binary representation of the rule.
    rule_b = np.array(
        [int(_) for _ in np.binary_repr(rule, 8)],
        dtype=np.int8)
    x = np.zeros((steps, size), dtype=np.int8)
    # Random initial state.
    x[0][int(len(x[0])/2)] = 1
    # Apply the step function iteratively.
    for i in range(steps - 1):
        x[i + 1, :] = step(x[i, :], rule_b)
    return x


fig = plt.figure(figsize=(8, 8))
rule = 30
x = generate(rule)
plt.imshow(x, interpolation='none',
          cmap=plt.cm.binary)
plt.title(str(rule))
