import click
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

@click.command()
@click.option('--h5-path', '-h', default='output/all_steps.h5', help='HDF5 file path')
@click.option('--mp4-path', '-m', default='nbody_h5_simulation.mp4', help='Output mp4 file path')
def main(h5_path, mp4_path):
    with h5py.File(h5_path, "r") as f:
        steps_group = f["/steps"]
        step_keys = sorted(steps_group.keys())
        num_steps = len(step_keys)
        positions = []
        for step in step_keys:
            pos = steps_group[step]["positions"][:]
            positions.append(pos)
        positions = np.array(positions)  # shape: (num_steps, N, 3)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    lim = 1.2e11
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    scat = ax.scatter([], [], [], s=20)
    ax.set_title("N-body Simulation (HDF5 -> MP4)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    def update(frame):
        pos = positions[frame]
        scat._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        ax.set_title(f"Step: {frame + 1}/{num_steps}")
        return scat

    ani = FuncAnimation(fig, update, frames=num_steps, interval=30, blit=False)
    ani.save(mp4_path, writer="ffmpeg", fps=30)

if __name__ == "__main__":
    main()