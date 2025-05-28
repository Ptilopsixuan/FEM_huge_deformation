from mylib import type, reader, writer
import numpy as np
import os

model_name = "beam"
input_path = os.path.join(os.getcwd(), "input", f"{model_name}.txt")
model = reader.readFile(input_path)
output = model.iterate()
output_path = os.path.join(os.getcwd(), "output", f"{model_name}.txt")
if not os.path.exists(os.path.dirname(output_path)):
    os.makedirs(os.path.dirname(output_path))
for i, o in enumerate(output):
    writer.writeFile(output_path.replace(".txt", f"_{i:03d}.txt"), o)





