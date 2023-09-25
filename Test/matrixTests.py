import numpy as np


numRow = 5
numCol = numRow


numMatrices = 10000

for i in range(numMatrices):
    # A = (np.random.randint(20, size=(numRow, numCol))) - 10
    A = np.random.uniform(-10, 10, size=(numRow, numCol))
    Ainv = np.linalg.inv(A)
    if i < 1:
        B = np.concatenate((A, Ainv), axis=1)
    else:
        B = np.concatenate((B, np.concatenate((A, Ainv), axis=1)), axis=0)

np.savetxt("Test/test.csv", B, delimiter=",", header=str(numRow))
