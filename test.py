from mylib import create, type, reader, writer
import numpy as np
import os

x = None  # Example x-coordinates
y = None  # Example y-coordinates
# theta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Example rotation angles

steps = 20000
for i in range(steps):
    file_str = create.create_beam(11, x, y, steps)
    input_path = os.path.join(os.getcwd(), "input","iterator", f"beam_{i:05d}.txt")
    if not os.path.exists(os.path.dirname(input_path)):
        os.makedirs(os.path.dirname(input_path))
    if i%10 == 0:
        print(f"Creating file: {input_path}")
    with open(input_path, "w") as f:
        f.write(file_str)

    model_name = f"beam_{i:05d}"  # Change this to the name of your model file without extension
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
    if i > 0.9 * steps and x[-1] >= 0:
        print(f"Beam {i} has exceeded the limit, stopping iteration.")
        break

# result = [float(r) for r in data]
# print(result)
# print(model.Pe)ste