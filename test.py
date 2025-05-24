from mylib import type, reader, writer

model_name = "3"
model = reader.readFile(f"input\\{model_name}.txt", shear = False)
model.integrateKe()
model.integratePe()
data = model.calculateA()
model.calculatePe()
model.transback()

output = type.OutputData(model.points, model.units, model.Pe)
writer.writeFile(f"output\\{model_name}.txt", output)

# result = [float(r) for r in data]
# print(result)
# print(model.Pe)