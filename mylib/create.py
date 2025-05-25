import numpy as np
import os

def create_beam(point_num: int = 11, x:list[float] = None, y:list[float] = None, steps: int = 1000) -> str:
    """
    Create a beam model with specified parameters.
    Args:
        x (list[float], optional): List of x-coordinates for the beam points. Defaults to None.
        y (list[float], optional): List of y-coordinates for the beam points. Defaults to None.
    Returns:
        str: A string representation of the beam model data.
    """
    E, v = 2e11, 0.3  # Young's modulus and Poisson's ratio
    b, h = 0.2, 0.3  # Width and height of the beam
    shape = 1 # Shape of the beam (1 for rectangular, 0 for circular)
    point_num = point_num  # Number of points along the beam
    if point_num < 2:
        raise ValueError("point_num must be at least 2.")
    length = 10 # Length of the beam
    goal_M = 2*np.pi*E*b*h **3/12/length  # Moment of inertia for rectangular section
    M = goal_M / steps  # Moment applied at the end of the beam
    # print(f"goal_M: {goal_M}")

    unit_num = point_num - 1
    # point_x = [i * length / unit_num for i in range(point_num)]
    # point_y = [0] * point_num
    if x is None or y is None:
        point_x = [i * length / unit_num for i in range(point_num)]
        point_y = [0] * point_num
    else:
        if len(x) != point_num or len(y) != point_num:
            raise ValueError("x and y must have the same length as point_num.")
        point_x = x
        point_y = y
    point_str = f"{point_num} #Points\n"
    for i in range(point_num):
        point_str += f"{i + 1} {point_x[i]} {point_y[i]}\n"

    material_str = "1 #Material\n"
    material_str += f"1 {E} {v}\n"

    plane_str = "1 #Plane\n"
    plane_str += f"1 {shape} {b} {h}\n"

    unit_str = f"{unit_num} #Unit\n"
    for i in range(unit_num):
        unit_str += f"{i + 1} {i + 1} {i + 2} 1 1\n"

    load_str = "1 #Load\n"
    load_str += f"{point_num} 0 0 {M}\n"

    constraint_str = "3 #Constraint\n"
    constraint_str += "1 1 0\n1 2 0\n1 3 0\n"

    hitch_str = "0 #Hitch\n"

    file_str = point_str + material_str + plane_str + unit_str + load_str + constraint_str + hitch_str
    return file_str
# print("point_x:", point_x)
# print("point_y:", point_y)

if __name__ == "__main__":
    file_str = create_beam()
    file_path = os.path.join(os.getcwd(), "input", "iterator", "beam.txt")
    with open(file_path, "w") as f:
        f.write(file_str)
    print(f"Beam data created and saved to {file_path}")
