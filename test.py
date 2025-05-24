from mylib import create, type, reader, writer
import numpy as np
import os

x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Example x-coordinates
y = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Example y-coordinates
# theta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Example rotation angles

steps = 2000
for i in range(steps):
    file_str = create.create_beam(x, y, steps)
    input_path = os.path.join(os.getcwd(), "input","iterator", f"beam_{i}.txt")
    if not os.path.exists(os.path.dirname(input_path)):
        os.makedirs(os.path.dirname(input_path))
    print(f"Creating file: {input_path}")
    with open(input_path, "w") as f:
        f.write(file_str)

    model_name = f"beam_{i}"  # Change this to the name of your model file without extension
    model = reader.readFile(input_path, shear = False)
    model.integrateKe()
    model.integratePe()
    data = model.calculateA()
    model.calculatePe()
    model.transback()

    output = type.OutputData(model.points, model.units, model.Pe)
    x = [point.x + point.ax for point in model.points]
    y = [point.y + point.ay for point in model.points]
    # theta = [point.theta + point.theta for point in model.points]

    output_path = os.path.join(os.getcwd(), "output","iterator", f"{model_name}.txt")
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    writer.writeFile(output_path, output)

    # in case of over
    if x[-1] >= 0:
        print(f"Beam {i} has exceeded the limit, stopping iteration.")
        break

# result = [float(r) for r in data]
# print(result)
# print(model.Pe)