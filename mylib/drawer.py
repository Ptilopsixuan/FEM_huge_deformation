from . import type
import matplotlib.pyplot as plt
import numpy as np

def d(data, path):
    fig, ax = plt.subplots(figsize=(8, 10), tight_layout = True)
    ax.set_aspect('equal')
    
    x, y = [],[]
    for p in data.points:
        x.append(p.x + p.ax)
        y.append(p.y + p.ay)
    scatter = plt.scatter(
    x,  # X坐标
    y,  # Y坐标
    c='blue',           # 点颜色
    s=40,               # 点大小
    alpha=0.7,          # 透明度
    edgecolors='black', # 边框颜色
    marker='o'          # 点形状
    )
    x_min = min(x) - 0.2
    y_min = min(y) - 0.2
    x_max = max(x) + 0.2
    y_max = max(y) + 0.2
    ax.set_xlim(x_min,x_max)
    ax.set_ylim(y_min,y_max)
    # plt.show()
    fig.savefig(path)
    plt.close()