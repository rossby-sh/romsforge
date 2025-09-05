from utils import ztosigma_numba
import numpy as np
import time

# 임시 테스트용 데이터
z = np.random.rand(30, 100, 100)
var = np.random.rand(30, 100, 100)
depth = np.linspace(-5000, 0, 40)

# 첫 호출 (컴파일)
t0 = time.time()
ztosigma_numba(var, z, depth)
t1 = time.time()

# 두 번째 호출 (캐시 사용)
ztosigma_numba(var, z, depth)
t2 = time.time()

print("First call (JIT compile):", t1 - t0)
print("Second call (cached):", t2 - t1)
