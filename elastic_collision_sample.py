import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

r1, r2 = 0.1, 0.1
pos1 = np.array([0.0, 0.0])
pos2 = np.array([1.5, 0.0])
vel1 = np.array([1.0, 0.1])   # Mavi cisme küçük y bileşeni
vel2 = np.array([-0.5, 0.0])
m1 = 1.5
m2 = 1.0

dt = 0.005


def elastic_collision_2d(m1, m2, v1, v2, x1, x2):
    delta_pos = x2 - x1
    dist = np.linalg.norm(delta_pos)
    n_hat = delta_pos / dist

    v1_n = np.dot(v1, n_hat)
    v2_n = np.dot(v2, n_hat)

    # Elastik çarpışma hız formülleri
    v1_n_new = (v1_n * (m1 - m2) + 2 * m2 * v2_n) / (m1 + m2)
    v2_n_new = (v2_n * (m2 - m1) + 2 * m1 * v1_n) / (m1 + m2)

    # Hız vektörlerini güncelle
    v1_new = v1 + (v1_n_new - v1_n) * n_hat
    v2_new = v2 + (v2_n_new - v2_n) * n_hat

    return v1_new, v2_new


fig, ax = plt.subplots()
ax.set_xlim(-1, 2.5)
ax.set_ylim(-1, 1)
ax.set_aspect('equal')
ball1 = plt.Circle(pos1, r1, fc='blue')
ball2 = plt.Circle(pos2, r2, fc='red')
ax.add_patch(ball1)
ax.add_patch(ball2)


def update(frame):
    global pos1, pos2, vel1, vel2

    pos1 += vel1 * dt
    pos2 += vel2 * dt

    dist = np.linalg.norm(pos2 - pos1)
    if dist <= r1 + r2:
        overlap = r1 + r2 - dist
        direction = (pos2 - pos1) / dist
        pos1 -= direction * overlap / 2
        pos2 += direction * overlap / 2

        vel1, vel2 = elastic_collision_2d(m1, m2, vel1, vel2, pos1, pos2)

        print(f"Çarpışma! frame={frame}")
        print(f"vel1 (mavi): {vel1}")
        print(f"vel2 (kırmızı): {vel2}")

    ball1.center = pos1
    ball2.center = pos2
    return ball1, ball2


ani = FuncAnimation(fig, update, frames=800, interval=10, blit=True)
plt.title("2D Tam Elastik Çarpışma Simülasyonu")
plt.show()
